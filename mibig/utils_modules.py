from mibig.errors import MibigError

class aSModule:
    def __init__(self, domains: list[str] = None, \
                 locus_tag: str = None, \
                 starterModule: bool | None = None, \
                 final_module: bool | None = None, \
                 locations: list[int] | None = None, \
                 ):
        if not all([domains, locus_tag]):
            raise MibigError("Domain list and locus tag are required")
        if not (location is None or len(location) == 2):
            raise MibigError("Location should be a list of two integers")

        self.domains = domains
        self.locus_tag = locus_tag
        self.starterModule = starterModule
        self.final_module = final_module
        self.location = location
