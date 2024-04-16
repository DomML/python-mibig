from typing import Any, Self

from mibig.converters.shared.common import Citation, validate_citation_list
from mibig.errors import ValidationError, ValidationErrorInfo


class FunctionEvidence:
    method: str
    references: list[Citation]

    VALID_METHODS = (
        "Other in vivo study",
        "Heterologous expression",
        "Knock-out",
        "Activity assay"
    )

    def __init__(self, method: str, references: list[Citation], validate: bool = True):
        self.method = method
        self.references = references

        if not validate:
            return

        errors = self.validate()
        if errors:
            raise ValidationError(errors)

    def validate(self) -> list[ValidationErrorInfo]:
        errors = []
        if self.method not in self.VALID_METHODS:
            errors.append(ValidationErrorInfo("method", f"Invalid method: {self.method}"))
        errors.extend(validate_citation_list(self.references))
        return errors

    @classmethod
    def from_json(cls, raw: dict[str, Any], **kwargs) -> Self:
        return cls(
            method=raw["method"],
            references=[Citation.from_json(ref) for ref in raw["references"]],
            **kwargs,
        )

    def to_json(self) -> dict[str, Any]:
        return {"method": self.method, "references": [ref.to_json() for ref in self.references]}


class MutationPhenotype:
    phenotype: str
    details: str | None
    references: list[Citation]

    def __init__(self, phenotype: str, references: list[Citation], details: str | None = None, validate: bool = True, **kwargs):
        self.phenotype = phenotype
        self.details = details
        self.references = references

        if not validate:
            return

        errors = self.validate(**kwargs)
        if errors:
            raise ValidationError(errors)

    def validate(self, **kwargs) -> list[ValidationErrorInfo]:
        errors = []
        if not self.phenotype:
            errors.append(ValidationErrorInfo("MutationPhenotype.phenotype", "Phenotype must be provided"))
        errors.extend(validate_citation_list(self.references))
        return errors

    @classmethod
    def from_json(cls, raw: dict[str, Any], **kwargs) -> Self:
        return cls(
            phenotype=raw["phenotype"],
            details=raw.get("details"),
            references=[Citation.from_json(ref) for ref in raw["references"]],
            **kwargs,
        )

    def to_json(self) -> dict[str, Any]:
        ret = {
            "phenotype": self.phenotype,
            "references": [ref.to_json() for ref in self.references],
        }
        if self.details:
            ret["details"] = self.details

        return ret


class GeneFunction:
    function: str
    evidence: list[FunctionEvidence]
    details: str | None
    mutation_phenotype: MutationPhenotype | None

    VALID_FUNCTIONS = (
        "Activation / processing",
        "Maturation",
        "Precursor",
        "Precursor biosynthesis",
        "Regulation",
        "Resistance/immunity",
        "Scaffold biosynthesis",
        "Tailoring",
        "Transport",
        "Other",
    )

    def __init__(self, function: str, evidence: list[FunctionEvidence],
                 details: str | None = None,
                 mutation_phenotype: MutationPhenotype | None = None,
                 validate: bool = True,
                 **kwargs
                ):
        self.function = function
        self.evidence = evidence
        self.details = details
        self.mutation_phenotype = mutation_phenotype

        if not validate:
            return

        errors = self.validate(**kwargs)
        if errors:
            raise ValidationError(errors)

    def validate(self, **kwargs) -> list[ValidationErrorInfo]:
        errors = []
        if self.function not in self.VALID_FUNCTIONS:
            errors.append(ValidationErrorInfo("GeneFunction.function", f"Invalid function: {self.function}"))

        if self.function == "Other" and not self.details:
            errors.append(ValidationErrorInfo("GeneFunction.details", "Details must be provided for 'Other' function"))

        if self.mutation_phenotype:
            errors.extend(self.mutation_phenotype.validate(**kwargs))

        for evidence in self.evidence:
            errors.extend(evidence.validate(**kwargs))
        return errors

    @classmethod
    def from_json(cls, raw: dict[str, Any], **kwargs) -> Self:
        return cls(
            function=raw["function"]["name"],
            evidence=[FunctionEvidence.from_json(evidence) for evidence in raw["evidence"]],
            details=raw["function"].get("details"),
            mutation_phenotype=MutationPhenotype.from_json(raw["mutation_phenotype"]) if "mutation_phenotype" in raw else None,
            **kwargs,
        )

    def to_json(self) -> dict[str, Any]:
        function = {
            "name": self.function,
        }
        if self.details:
            function["details"] = self.details

        ret = {
            "function": function,
            "evidence": [evidence.to_json() for evidence in self.evidence],
        }
        if self.mutation_phenotype:
            ret["mutation_phenotype"] = self.mutation_phenotype.to_json()

        return ret

