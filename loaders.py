import pandas as pd


class PubMedLoader:
    def load(
        self,
        pubmed_df: pd.DataFrame,
        fmc_db_df: pd.DataFrame,
    ) -> pd.DataFrame:
        df = pubmed_df.copy()
        df["isStored"] = (df["firstName"] + " " + df["lastName"]).isin(
            fmc_db_df["lastName"].str.split(",").str.get(0)
        )
        df.to_json("output_file.json", orient="records", indent=4)

        print(df)
        return df
