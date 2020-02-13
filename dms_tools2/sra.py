"""
===================
sra
===================

Functions for downloading / handling data from the 
`Sequence Read Archive <https://www.ncbi.nlm.nih.gov/sra>`_ (SRA).
"""


import distutils.version
import multiprocessing
import os
import shutil
import subprocess


def fastqFromSRA(samples, fastq_dump, fastqdir, aspera=None,
        overwrite=False, passonly=True, no_downloads=False,
        ncpus=1):
    """Download data from SRA and extract FASTQ files.

    Currently only works for runs containing paired-end reads.

    Args:
        `samples` (`pandas.DataFrame`)
            A dataframe that must have columns named `run`
            and `name`. The `run` column gives SRA run accessions
            (e.g., `SRR5241726`), and the `name` column gives the
            name for the run used in the final FASTQ files.
            Will be modified to include `R1` and `R2` columns.
        `fastq_dump` (str)
            Path to ``fastq-dump`` executable. Requires a
            version >= 2.8.
        `fastqdir` (str)
            Directory in which to place the FASTQ files.
            Created if it does not already exist.
        `aspera` (`None` or 2-tuple)
            If `None`, use ``fastq-dump`` for downloads (this is slower).
            However, downloads are faster with 
            `aspera <https://www.ncbi.nlm.nih.gov/books/NBK242625/>`_
            To use aspera, specify the 2-tuple `(ascp, asperakey)` where 
            `ascp` is path to ``ascp`` executable, and `asperakey` is the
            `key <https://www.ncbi.nlm.nih.gov/sra/docs/aspera-key-pairs/>`_.
        `overwrite` (bool)
            If file already exists, do we overwrite it or just
            use the existing one? If `False` and all output files
            already exist, then nothing is done and `fastq_dump`
            no longer even needs to be a valid path.
        `passonly` (bool)
            Keep only reads with a passing `READ_FILTER` value.
        `no_downloads` (bool)
            If `True`, do **not** actually download the files,
            but instead just add the columns to `samples`.
        `ncpus` (int)
            Use this many CPUs to parallelize downloads, downgrades
            number if it exceeds max available.

    Result:
        Upon completion, the directory `fastqdir` contains
        files of the form ``<name>_R1.fastq.gz`` and
        ``<name>_R2.fastq.gz`` for all names in `samples`.
        These names have been added as the columns `R1` and
        `R2` to `samples`. Note that the file names but
        **not** the directory names are added.
    """
    assert {'run', 'name'} <= set(samples.columns), \
            '`samples` does not have columns `run` and `name`'
    samples['R1'] = samples['name'] + '_R1.fastq.gz'
    samples['R2'] = samples['name'] + '_R2.fastq.gz'

    if (all([os.path.isfile(os.path.join(fastqdir, f))
            for r in ['R1', 'R2'] for f in samples[r]])
            and not overwrite):
        return # all files already present, nothing to do

    if no_downloads:
        return # don't do anything

    assert shutil.which(fastq_dump), ("fastq-dump not installed in a "
            "location accessible with command {0}".format(fastq_dump))
    fastq_dump_version = (subprocess.check_output([fastq_dump, '--version'])
                          .decode('utf-8').split(':'))
    if len(fastq_dump_version) == 2:
        fastq_dump_version = fastq_dump_version[1].strip()
    else:
        fastq_dump_version = (subprocess.check_output([fastq_dump, '--help'])
                              .decode('utf-8').split(':')[-1].strip())
    fastq_dump_minversion = '2.8'
    if not (distutils.version.LooseVersion(fastq_dump_version) >=
            distutils.version.LooseVersion(fastq_dump_minversion)):
        raise RuntimeError("fastq-dump version {0} is installed. You need "
            "at least version {1}".format(fastq_dump_version, 
            fastq_dump_minversion))

    if aspera is not None:
        (ascp, asperakey) = aspera
        assert shutil.which(ascp), "ascp is not installed at {0}".format(ascp)
        assert os.path.isfile(asperakey), "no aspera key {0}".format(asperakey)

    if not os.path.isdir(fastqdir):
        os.mkdir(fastqdir)

    args_list = []
    expect_files = []
    for row in samples.iterrows():
        row = row[1]
        (run, name, r1base, r2base) = (
                row['run'], row['name'], row['R1'], row['R2'])
        r1 = os.path.join(fastqdir, r1base)
        r2 = os.path.join(fastqdir, r2base)
        expect_files += [r1, r2]

        if (not overwrite) and os.path.isfile(r1) and os.path.isfile(r2):
            continue # files already exist
        else:
            args_list.append((fastqdir,
                                run,
                                passonly,
                                aspera,
                                overwrite,
                                fastq_dump,
                                r1,
                                r2,
                                ))

    ncpus = min(ncpus, multiprocessing.cpu_count())
    if ncpus == 1:
        for args in args_list:
            _getSingleFastqFromSRA(*args)
    else:
        with multiprocessing.Pool(ncpus) as pool:
            pool.starmap(_getSingleFastqFromSRA, args_list)

    missing_files = [f for f in expect_files if not os.path.isfile(f)]
    if missing_files:
        raise RuntimeError('Failed to create these files:\n  ' +
                           '\n  '.join(missing_files))


def _getSingleFastqFromSRA(fastqdir,
                           run,
                           passonly,
                           aspera,
                           overwrite,
                           fastq_dump,
                           r1,
                           r2,
                           ):
    """Helper function for :func:`fastqFromSRA`."""
    # names of files produced by fastq-dump
    r1init = os.path.join(fastqdir, '{0}{1}_1.fastq.gz'.format(
                    run, {True:'_pass', False:''}[bool(passonly)]))
    r2init = os.path.join(fastqdir, '{0}{1}_2.fastq.gz'.format(
                    run, {True:'_pass', False:''}[bool(passonly)]))

    if aspera is not None:
        # use aspera to download
        (ascp, asperakey) = aspera

        # directory name explanation here: 
        # https://www.ncbi.nlm.nih.gov/books/NBK158899/
        target = ('anonftp@ftp.ncbi.nlm.nih.gov:'
                  '/sra/sra-instant/reads/ByRun/sra/SRR/{0}/{1}'
                  '/{1}.sra'.format(run[ : 6], run))

        asperacmds = [ascp, 
                      '--overwrite', 
                      {True:'always', False:'never'}[bool(overwrite)],
                      '-T',
                      '-q', # quiet mode
                      '--policy', 'fair',
                      '-l500m',
                      '-i', asperakey,
                      target,
                      fastqdir]
        asperaout = subprocess.check_call(asperacmds)
        downloadedrun = os.path.join(fastqdir, run + '.sra')
        assert os.path.isfile(downloadedrun)

    else:
        downloadedrun = None

    # run fastq-dump
    cmds = [fastq_dump,
            '--outdir', fastqdir,
            '--gzip',
            '--readids',
            '--skip-technical',
            '--dumpbase',
            '--clip',
            '--split-files',
            ]
    if passonly:
        cmds += ['--read-filter', 'pass']
    if downloadedrun:
        cmds.append(downloadedrun)
    else:
        cmds.append(run)
    out = subprocess.check_output(cmds)

    os.rename(r1init, r1)
    os.rename(r2init, r2)

    if downloadedrun:
        os.remove(downloadedrun)
                

if __name__ == '__main__':
    import doctest
    doctest.testmod()
