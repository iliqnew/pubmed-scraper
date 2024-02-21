pubmed-scraper is suposed to get data from the [PubMed API](https://www.ncbi.nlm.nih.gov/home/develop/api/) using [biopython package](https://biopython.org/wiki/Documentation)

Steps of the development:

Initially I looked on the NCBI website if they allow data to be scraped but then I found out about biopython.
I began reading the docs for it's my first time I work with it. There were some back-and-forwards but after some docs I finally managed to make it extract data.
I separated the code by functionallities and then tought of an architecture.

There are four main classes:
- Scrapers
- Transformers (Go Optimus!)
- Loaders
- Pipelines

Scraper classes are ment to hold the scraping business logic
Transformers - make transformations to the scraped data
Loaders - load the transformed data (ui, api, file delivery, etc.)
Pipelines' goal is to hold the solution for the desired result on a higher level describing the combination of scrapers, transformers and loaders needed

These four classes inherit from abstract ones which lay the foundation of a Strategy design pattern which is here used to ensure scalability

Don't hesitate to educate me about any topic that relates. Thanks in advance
