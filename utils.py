import json

from typing import Union, Literal, Any


with open("person_profiles.json", "rb") as f:
    STORED_PROFILES = json.loads(f.read())


def clean_names(
    fore_name: Union[str, None], last_name: Union[str, None]
) -> dict[Literal["firstName", "lastName"], Union[str, None]]:
    # apply data cleaning logic here

    return {"firstName": fore_name, "lastName": last_name}


def clean_is_stored(author: dict[str, Any]) -> dict[Literal["isStored"], bool]:
    fn = author["firstName"]
    ln = author["lastName"]

    _first_name = fn + " " if fn else ""
    _last_name = ln if ln else ""

    is_stored = _first_name + _last_name in [
        profile["lastName"].split(",")[0]
        for profile in STORED_PROFILES
        if not profile["lastName"] is None
    ]

    return {"isStored": is_stored}
