import json
import numpy as np

from abc import ABC, abstractmethod

from utils import search_query, get_publications_by_ids, clean_names


class Scraper(ABC):
    author_t = {
        # "_id": {"$oid": "5c6ad92e1fba1942ccb101c2"},
        "firstName": None,
        # "middleName": None,
        "lastName": None,  # "Tjalling W de Vries, dr.",
        "affiliation": None,
        # "degrees": None,
        # "phone": "0031582863390",
        # "phoneExt": None,
        # "email": "tjalling.de.vries@znb.nl",
        # "roleSpecified": True,
        # "role": 0,
        # "roles": [1],
        # "personType": 0,
        # "name": "Tjalling W de Vries, dr.",
        # "formattedName": "Tjalling De Vries",
        # "isMarkedAsDuplicate": False,
        # "lastDuplicationCheck": {"$date": "2022-04-26T01:44:11.340Z"},
        # "personTypeFixDate": {"$date": {"$numberLong": "-62135596800000"}},
        # "createdOn": {"$date": {"$numberLong": "-62135596800000"}},
        # "debarmentData": None,
        # "formattedAddress": "Leeuwarden, Friesland, Netherlands",
        # "website": None,
        # "placeId": None,
        # "location": None,
        # "contactsDownloadedFrom": None,
        # "googleSERPProcessed": False,
        # "notRealPerson": False,
        # "npi": None,
        # "medicalSchool": None,
        # "specialties": None,
        # "affiliateHospitals": None,
        # "medicareProcessed": True,
        # "medicareProcessedFix": True,
        # "phonesInfo": [
        #     {"trialId": "NCT03654508", "value": "0031582863390"},
        #     {"trialId": "NCT02197780", "value": "0031 582 866 666"},
        # ],
        # "emailsInfo": [
        #     {"trialId": "NCT03654508", "value": "tjalling.de.vries@znb.nl"},
        #     {"trialId": "NCT02197780", "value": "Tjalling.de.Vries@ZNB.nl"},
        # ],
        # "emailsInfoUpdated": True,
        # "fDAInspections": None,
        # "normalizedEmails": ["TJALLING.DE.VRIES@ZNB.NL"],
        # "normalizedPhones": ["0031582863390", "0031582866666"],
        # "trials": ["NCT03654508", "NCT02197780"],
    }

    def __init__(
        self,
        batch_size: int = 9999,
        limit: int = np.inf,
    ) -> None:
        self.batch_size = batch_size
        self.limit = limit

    @abstractmethod
    def scrape(self): ...

    @abstractmethod
    def parse(self): ...


class PubMedScraper(Scraper):
    query = "Ulcerative Colitis"

    def scrape(self):
        retstart = 0
        retmax = self.batch_size
        keep_fetching = True
        batch = 1

        while keep_fetching:
            studies_ids = search_query(self.query, retstart=retstart, retmax=retmax)
            studies_id_list = studies_ids["IdList"]

            # scrape the authors of n publications; self.limit = n
            if batch * self.batch_size == self.limit:
                keep_fetching = False
            elif batch * self.batch_size > self.limit:
                studies_id_list = studies_id_list[: self.limit % self.batch_size]
                keep_fetching = False

            raw_studies = get_publications_by_ids(
                studies_id_list,
            )
            for author in self.parse(raw_studies["PubmedArticle"]):
                yield author

            retstart = retstart + retmax
            batch += 1

    def parse(self, raw_input):
        for raw_study in raw_input:
            raw_authors = (
                raw_study["MedlineCitation"]["Article"].get("AuthorList") or []
            )
            for raw_author in raw_authors:
                author = self.author_t.copy()

                cleaned_names = clean_names(
                    raw_author.get("ForeName"), raw_author.get("LastName")
                )
                first_name = cleaned_names["firstName"]
                last_name = cleaned_names["lastName"]
                if first_name is last_name is None:
                    continue

                author.update({"firstName": first_name, "lastName": last_name})
                yield author


class FindMeCureDBScraper(Scraper):

    def scrape(self):
        with open("person_profiles.json", "rb") as f:
            profiles = json.loads(f.read())

        return profiles

    def parse(self, raw_input):
        pass
