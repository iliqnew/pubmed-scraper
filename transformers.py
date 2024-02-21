import pandas as pd
from abc import ABC, abstractmethod


class Transformer(ABC):
    @abstractmethod
    def transform(self, raw_input): ...


class PubMedTransformer(Transformer):

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        return df


class FindMeCureDBTransformer(Transformer):

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        return df
