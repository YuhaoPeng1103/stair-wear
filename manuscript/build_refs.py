"""Build refs.bib from verified DOIs via Crossref content negotiation.

Every DOI in DOI_LIST was resolved during preparation. The script prints the
retrieved title so a mismatch is visible. Book/standard entries are added
manually at the end.
"""
from __future__ import annotations

import time
import urllib.request
from pathlib import Path

OUT = Path(__file__).with_name("refs.bib")

DOI_LIST = [
    "10.1063/1.1721448",                 # Archard 1953
    "10.3390/app16126025",               # Yu et al. 2026
    "10.54097/yawxst60",                 # Li et al. 2025
    "10.26855/jamc.2025.12.006",         # Cao et al. 2025
    "10.1145/3745533.3745639",           # Lu et al. 2025
    "10.3233/ATDE250215",                # Cai et al. 2025
    "10.5220/0015109700005244",          # Mao et al. 2026
    "10.1093/bioinformatics/btp358",     # Raue et al. 2009
    "10.1111/febs.12276",                # Kreutz et al. 2013
    "10.1016/j.coisb.2021.03.005",       # Wieland et al. 2021
    "10.1016/S0140-6736(86)90837-8",     # Bland & Altman 1986
    "10.3390/rs3061104",                 # Remondino 2011
    "10.1016/S0966-6362(01)00162-X",     # Riener et al. 2002
    "10.1007/BF01068419",                # Schuirmann 1987
    "10.1016/j.culher.2014.06.007",      # Delgado Rodrigues 2015
    "10.1179/019713600806113176",        # Charola 2000
    "10.6028/NIST.TN.1839",              # Kuligowski et al. 2014
    "10.1016/j.trc.2023.104468",         # Pouw et al. 2024
]

MANUAL = r"""
@book{Johnson1985,
  author    = {Johnson, K. L.},
  title     = {Contact Mechanics},
  publisher = {Cambridge University Press},
  address   = {Cambridge},
  year      = {1985}
}

@book{Rabinowicz1995,
  author    = {Rabinowicz, Ernest},
  title     = {Friction and Wear of Materials},
  edition   = {2},
  publisher = {Wiley},
  address   = {New York},
  year      = {1995}
}

@book{Bhushan2013,
  author    = {Bhushan, Bharat},
  title     = {Principles and Applications of Tribology},
  edition   = {2},
  publisher = {Wiley},
  address   = {Hoboken, NJ},
  year      = {2013}
}

@book{DoehnePrice2010,
  author    = {Doehne, Eric and Price, Clifford A.},
  title     = {Stone Conservation: An Overview of Current Research},
  edition   = {2},
  publisher = {Getty Conservation Institute},
  address   = {Los Angeles},
  year      = {2010}
}

@book{SiegesmundSnethlage2014,
  editor    = {Siegesmund, Siegfried and Snethlage, Rolf},
  title     = {Stone in Architecture: Properties, Durability},
  edition   = {5},
  publisher = {Springer},
  address   = {Berlin},
  year      = {2014}
}

@book{BurnhamAnderson2002,
  author    = {Burnham, Kenneth P. and Anderson, David R.},
  title     = {Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach},
  edition   = {2},
  publisher = {Springer},
  address   = {New York},
  year      = {2002}
}

@book{Gelman2013,
  author    = {Gelman, Andrew and Carlin, John B. and Stern, Hal S. and Dunson, David B. and Vehtari, Aki and Rubin, Donald B.},
  title     = {Bayesian Data Analysis},
  edition   = {3},
  publisher = {CRC Press},
  address   = {Boca Raton, FL},
  year      = {2013}
}

@book{EfronTibshirani1993,
  author    = {Efron, Bradley and Tibshirani, Robert J.},
  title     = {An Introduction to the Bootstrap},
  publisher = {Chapman and Hall/CRC},
  address   = {New York},
  year      = {1993}
}

@misc{ICOMOS2013,
  author       = {{ICOMOS}},
  title        = {The Burra Charter: The Australia ICOMOS Charter for Places of Cultural Significance},
  year         = {2013},
  howpublished = {Australia ICOMOS, Burwood}
}

@misc{HistoricEngland2017,
  author       = {{Historic England}},
  title        = {Photogrammetric Applications for Cultural Heritage: Guidance for Good Practice},
  year         = {2017},
  howpublished = {Historic England, Swindon},
  url          = {https://historicengland.org.uk/images-books/publications/photogrammetric-applications-for-cultural-heritage/}
}

@misc{ISO5725,
  author       = {{International Organization for Standardization}},
  title        = {ISO 5725-2:2019 Accuracy (Trueness and Precision) of Measurement Methods and Results --- Part 2},
  year         = {2019},
  howpublished = {ISO, Geneva}
}

@misc{Tribunnews2015,
  author       = {{Tribunnews}},
  title        = {Cegah keausan batu, pengunjung Borobudur diharuskan pakai sandal khusus},
  year         = {2015},
  howpublished = {13 October 2015},
  url          = {https://www.tribunnews.com/regional/2015/10/13/cegah-keausan-batu-pengunjung-borobudur-diharuskan-pakai-sandal-khusus},
  note         = {Grey literature; reports a Borobudur Conservation Office projection}
}

@misc{JCGM2008,
  author       = {{Joint Committee for Guides in Metrology}},
  title        = {JCGM 100:2008 Evaluation of Measurement Data --- Guide to the Expression of Uncertainty in Measurement},
  year         = {2008},
  howpublished = {BIPM, S\`evres}
}
"""


def fetch(doi: str) -> str | None:
    req = urllib.request.Request(
        f"https://doi.org/{doi}",
        headers={"Accept": "application/x-bibtex", "User-Agent": "ref-builder/1.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=25) as r:
            return r.read().decode("utf-8", errors="replace").strip()
    except Exception as exc:  # noqa: BLE001
        print(f"FAILED {doi}: {exc}")
        return None


def main() -> None:
    entries = []
    for doi in DOI_LIST:
        txt = fetch(doi)
        if txt:
            entries.append(txt)
            title = ""
            for line in txt.splitlines():
                if "title=" in line:
                    title = line.split("title=", 1)[1][:120]
                    break
            print(f"OK {doi}  {title}")
        time.sleep(0.2)
    OUT.write_text("\n\n".join(entries) + "\n" + MANUAL, encoding="utf-8")
    print(f"\nWrote {OUT} with {len(entries)} Crossref entries + manual books/standards.")


if __name__ == "__main__":
    main()
