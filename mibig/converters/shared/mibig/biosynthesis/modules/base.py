from enum import StrEnum
from typing import Any, Self

from mibig.converters.shared.common import (
    Citation,
    Evidence,
    GeneId,
    QualityLevel,
)
from mibig.converters.shared.mibig.biosynthesis.common import Monomer
from mibig.errors import ValidationError, ValidationErrorInfo

from .cal import CAL
from .core import ModuleInfo
from .nrps import NrpsTypeI
from .other import Other
from .pks import (
    PksIterative,
    PksModular,
    PksModularStarter,
    PksTransAt,
    PksTransAtStarter,
)

ExtraInfo = ModuleInfo


class NcaEvidence(Evidence):
    VALID_METHODS = (
        "Sequence-based prediction",
        "Structure-based inference",
        "Activity assay",
    )


class NonCanonicalActivity:
    evidence: list[NcaEvidence]
    iterations: int | None = None
    non_elongating: bool | None = None
    skipped: bool | None = None

    def __init__(
        self,
        evidence: list[NcaEvidence],
        iterations: int | None = None,
        non_elongating: bool | None = None,
        skipped: bool | None = None,
        validate: bool = True,
        **kwargs,
    ) -> None:
        self.evidence = evidence
        self.iterations = iterations
        self.non_elongating = non_elongating
        self.skipped = skipped

        if not validate:
            return

        errors = self.validate(**kwargs)
        if errors:
            raise ValidationError(errors)

    def validate(self, **kwargs) -> list[ValidationErrorInfo]:
        errors = []

        quality: QualityLevel | None = kwargs.get("quality")

        if quality != QualityLevel.QUESTIONABLE:
            if not self.evidence:
                errors.append(
                    ValidationErrorInfo(
                        "NonCanonicalActivity.evidence",
                        "At least one evidence must be provided",
                    )
                )
        for evidence in self.evidence:
            errors.extend(evidence.validate(**kwargs))

        return errors

    @classmethod
    def from_json(cls, raw: dict[str, Any], **kwargs) -> Self:
        return cls(
            evidence=[NcaEvidence.from_json(evidence, **kwargs) for evidence in raw["evidence"]],
            iterations=raw.get("iterations"),
            non_elongating=raw.get("nonElongating"),
            skipped=raw.get("skipped"),
            **kwargs,
        )

    def to_json(self) -> dict[str, Any]:
        ret: dict[str, Any] = {
            "evidence": [evidence.to_json() for evidence in self.evidence],
        }
        if self.iterations is not None:
            ret["iterations"] = self.iterations
        if self.non_elongating is not None:
            ret["nonElongating"] = self.non_elongating
        if self.skipped is not None:
            ret["skipped"] = self.skipped
        return ret


class ModuleType(StrEnum):
    CAL = "cal"
    NRPS_TYPE1 = "nrps-type1"
    NRPS_TYPE6 = "nrps-type6"
    OTHER = "other"
    PKS_ITERATIVE = "pks-iterative"
    PKS_MODULAR = "pks-modular"
    PKS_TRANS_AT = "pks-trans-at"
    PKS_MODULAR_STARTER = "pks-modular-starter"
    PKS_TRANS_AT_STARTER = "pks-trans-at-starter"


MAPPING = {
    ModuleType.CAL: CAL,
    ModuleType.NRPS_TYPE1: NrpsTypeI,
    ModuleType.NRPS_TYPE6: NrpsTypeI,  # There's basically no difference between the two
    ModuleType.OTHER: Other,
    ModuleType.PKS_ITERATIVE: PksIterative,
    ModuleType.PKS_MODULAR: PksModular,
    ModuleType.PKS_TRANS_AT: PksTransAt,
    ModuleType.PKS_MODULAR_STARTER: PksModularStarter,
    ModuleType.PKS_TRANS_AT_STARTER: PksTransAtStarter,
}


class Module:
    module_type: ModuleType
    name: str
    genes: list[GeneId]
    active: bool
    extra_info: ExtraInfo
    integrated_monomers: list[Monomer]
    non_canonical_activity: NonCanonicalActivity | None
    comment: str | None

    def __init__(
        self,
        module_type: ModuleType,
        name: str,
        genes: list[GeneId],
        active: bool,
        extra_info: ExtraInfo,
        integrated_monomers: list[Monomer],
        non_canonical_activity: NonCanonicalActivity | None = None,
        comment: str | None = None,
        validate: bool = True,
        **kwargs,
    ) -> None:
        self.module_type = module_type
        self.name = name
        self.genes = genes
        self.active = active
        self.extra_info = extra_info
        self.integrated_monomers = integrated_monomers
        self.non_canonical_activity = non_canonical_activity
        self.comment = comment

        if not validate:
            return

        errors = self.validate(**kwargs)
        if errors:
            raise ValidationError(errors)

    def validate(self, **kwargs) -> list[ValidationErrorInfo]:
        errors = []

        errors.extend(self.extra_info.validate(**kwargs))

        for monomer in self.integrated_monomers:
            errors.extend(monomer.validate())

        for gene in self.genes:
            errors.extend(gene.validate(record=kwargs.get("record")))

        if self.non_canonical_activity:
            errors.extend(self.non_canonical_activity.validate(**kwargs))

        return errors

    @property
    def references(self) -> list[Citation]:
        references = set()
        for monomer in self.integrated_monomers:
            references.update(monomer.references)
        for domain in self.extra_info.get_domains():
            references.update(domain.references)
        return sorted(list(references))

    @classmethod
    def from_json(cls, raw: dict[str, Any], **kwargs) -> Self:
        module_type = ModuleType(raw["type"])
        extra_info = MAPPING[ModuleType(module_type)].from_json(raw, **kwargs)
        nc_activity = None
        if raw.get("non_canonical_activity"):
            nc_activity = NonCanonicalActivity.from_json(raw["non_canonical_activity"], **kwargs)

        return cls(
            module_type=module_type,
            name=raw["name"],
            genes=[GeneId.from_json(gene, **kwargs) for gene in raw["genes"]],
            active=raw["active"],
            extra_info=extra_info,
            comment=raw.get("comment"),
            integrated_monomers=[
                Monomer.from_json(monomer)
                for monomer in raw.get("integrated_monomers", [])
            ],
            non_canonical_activity=nc_activity,
            **kwargs,
        )

    def to_json(self) -> dict[str, Any]:
        ret = {
            "type": self.module_type.value,
            "name": self.name,
            "genes": [gene.to_json() for gene in self.genes],
            "active": self.active,
            **self.extra_info.to_json(),
        }
        if self.integrated_monomers:
            ret["integrated_monomers"] = [
                monomer.to_json() for monomer in self.integrated_monomers
            ]

        if self.non_canonical_activity:
            ret["non_canonical_activity"] = self.non_canonical_activity.to_json()

        if self.comment:
            ret["comment"] = self.comment

        return ret
