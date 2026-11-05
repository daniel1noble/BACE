# Bird illustrations for the manuscript

Four transparent scientific illustrations for panel d of the manuscript's case-study figure. Generated with the built-in `image_gen` tool on 2026-11-05; exact prompts are in `generation-prompts.json`. The illustrations are illustrative, not measurements, and are not drawn to scale.

The examples were selected from **observed** records in `../../bird_social/data/social_data.rds`. Each selected species has exactly one positive social-grouping indicator and no lek flag. The category is the largest reported grouping under this dataset's coding; it does not imply that the species never exhibits other behaviours.

| Category | Species | Asset |
|---|---|---|
| Solitary | Common kingfisher, *Alcedo atthis* | `common_kingfisher.png` |
| Singly and pairs | Western barn owl, *Tyto alba* | `western_barn_owl.png` |
| Pairs and family groups | Laughing kookaburra, *Dacelo novaeguineae* | `laughing_kookaburra.png` |
| Flocks or colonies | Sulphur-crested cockatoo, *Cacatua galerita* | `sulphur_crested_cockatoo.png` |

Appearance and nomenclature references:

- [Common kingfisher: Naturalis](https://birds-europe.linnaeus.naturalis.nl/linnaeus_ng/app/views/species/taxon.php?id=130758)
- [Barn owl: Natural History Museum](https://www.nhm.ac.uk/discover/barn-owl-tyto-alba.html); [Western barn owl: Birds of the World](https://doi.org/10.2173/bow.webowl1.01)
- [Laughing kookaburra: Australian Museum](https://australian.museum/learn/animals/birds/laughing-kookaburra/)
- [Sulphur-crested cockatoo: Australian Museum](https://australian.museum/learn/animals/birds/sulphur-crested-cockatoo/)

All four PNGs retain their original RGBA transparency and were copied byte-for-byte from the generated assets. The plotting code in `ms/ms.qmd` places them in the category-label margin, outside the data coordinates, with a fixed square aspect ratio. The category-to-species mapping was checked against the saved input data, and the final figure was visually inspected.
