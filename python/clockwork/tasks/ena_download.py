from clockwork import ena_downloader

def run(options):
    downloader = ena_downloader.EnaDownloader(
        options.data_tsv,
        options.output_dir,
        options.site,
        options.lab,
        options.date,
        options.dataset_name,
        download_threads=options.download_threads,
        md5_threads=options.md5_threads,
    )
    downloader.run()

