import inspect
from typing import Any, Dict
import attr
from daesim.management import ManagementModule  # or PlantGrowthPhases, etc.

class DAESIMModulePartialLoad:
    """
    Partial-like callable wrapper for DAESIM modules: captures default constructor
    arguments (from attrs or signature) and allows overrides when called.
    """
    def __init__(self, module_cls: Any, **user_args):
        # Determine target class
        self._module_cls = module_cls

        # Attempt to pull defaults from an attrs-based instance
        defaults: Dict[str, Any] = {}
        try:
            instance = module_cls() if isinstance(module_cls, type) else module_cls
            defaults = attr.asdict(instance)
        except Exception:
            # Fallback to signature defaults
            sig = inspect.signature(self._module_cls.__init__)
            for name, param in sig.parameters.items():
                if name == 'self':
                    continue
                if param.default is not inspect._empty:
                    defaults[name] = param.default

        defaults.update(user_args)
        self._partial_kwargs = defaults

    def __call__(self, **user_args) -> Any:
        kwargs = self._partial_kwargs.copy()
        kwargs.update(user_args)
        sig = inspect.signature(self._module_cls.__init__)
        valid = {name for name in sig.parameters if name != 'self'}
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in valid}
        return self._module_cls(**filtered_kwargs)

    def to_dict(self) -> Dict[str, Any]:
        """
        Return the captured default kwargs (before call-time overrides).
        """
        return dict(self._partial_kwargs)