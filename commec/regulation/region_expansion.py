"""
Temporary storage for which countries are in the EU vs the AG.
In the future we will likely want to import this info from some other dynamic source.
"""
from commec.regulation.containers import Region

region_codes = {
    (Region("EuropeanUnion", "EU"), ["AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", 
                                     "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", 
                                     "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO", 
                                     "SE", "SI", "SK"]),
    (Region("AustraliaGroup", "AG"), ["AR", "AT", "AU", "BE", "BG", "CA", "CH", "CY", 
                                      "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GB", 
                                      "GR", "HR", "HU", "IE", "IN", "IS", "IT", "JP", 
                                      "KR", "LT", "LU", "LV", "MT", "MX", "NL", "NO", 
                                      "NZ", "PL", "PT", "RO", "SE", "SI", "SK", "TR", 
                                      "UA", "US"])
}