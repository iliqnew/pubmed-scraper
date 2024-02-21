import pandas as pd

from abc import ABC, abstractmethod

from scrapers import PubMedScraper, FindMeCureDBScraper

# from transformers import PubMedTransformer, FindMeCureDBTransformer

from loaders import PubMedLoader


class Pipeline(ABC):
    @abstractmethod
    def run(): ...


class PubMedPipeline(Pipeline):
    def run():
        pubmed_results = PubMedScraper(
            batch_size=2500,  # max 9999 publications
            limit=5000,
        ).scrape()
        fmc_db_results = FindMeCureDBScraper().scrape()

        pubmed_df = pd.json_normalize(pubmed_results)
        fmc_db_df = pd.json_normalize(fmc_db_results)

        # pubmed_df = PubMedTransformer().transform(pubmed_df)
        # fmc_db_df = FindMeCureDBTransformer().transform(fmc_db_df)

        return PubMedLoader().load(pubmed_df, fmc_db_df)
