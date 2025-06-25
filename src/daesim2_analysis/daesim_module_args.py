

import inspect
from typing import Any, Dict
import attr

class DAESIMModuleArgs:
    """
    Wrapper to merge default attrs-based or plain-object module args with user overrides,
    then instantiate only with valid constructor parameters.
    """
    def __init__(self, module: Any, **user_args):
        self.module = module
        if isinstance(module, type):
            module = module()

        # Extract defaults from attrs or __dict__
        try:
            source = attr.asdict(module)
        except (attr.exceptions.NotAnAttrsClassError, TypeError, ValueError):
            source = vars(module)

        # Build default_args dict: prefer user_args over module defaults, skip None
        self.default_args: Dict[str, Any] = {}
        for name, value in source.items():
            if not name.startswith('_') and value is not None:
                self.default_args[name] = user_args.get(name, value)
        # Allow extra user_args not in defaults
        for name, value in user_args.items():
            if name not in self.default_args:
                self.default_args[name] = value

    def to_dict(self) -> Dict[str, Any]:
        """
        Return the merged default and user override arguments.
        """
        return dict(self.default_args)

    def create(self, **user_kwargs) -> Any:
        """
        Instantiate the module, merging defaults, overrides, and filtering to valid __init__ params.
        """
        # Merge existing defaults with new overrides
        merged = self.default_args.copy()
        merged.update(user_kwargs)

        # Determine class to instantiate
        cls = self.module if isinstance(self.module, type) else type(self.module)

        # Filter to valid constructor parameters
        sig = inspect.signature(cls.__init__)
        valid = {p for p in sig.parameters if p != 'self'}
        filtered = {k: v for k, v in merged.items() if k in valid}

        return cls(**filtered)
