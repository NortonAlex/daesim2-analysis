from typing import Any, Dict
from daesim.management import ManagementModule
from pprint import pprint
import attr

class DAESIMModuleArgs:
    def __init__(self, module):
        self.module = module
        if isinstance(module, type):
            module = module()
        try:
            source = attr.asdict(module)
        except (attr.exceptions.NotAnAttrsClassError, TypeError, ValueError):

            source = vars(module)

        public_items = {name: value for name, value in source.items() if not name.startswith('_')}
        for name, value in public_items.items():
            if value is not None:
                setattr(self, name, value)

    def to_dict(self) -> Dict[str, Any]:
        """
        Return a dictionary of all public management arguments.
        """
        return {name: getattr(self, name) for name in vars(self) if not name.startswith('_')}

    def merge_current_and_default_args(self, user_args: dict):
        current_args = self.to_dict()
        merged_args = {}
        for arg in current_args:
            if arg in user_args:
                merged_args[arg] = user_args[arg]
            else:
                merged_args[arg] = current_args[arg]
        return merged_args

    def create(self, **kwargs):
        args = self.merge_current_and_default_args(kwargs)
        module = self.module
        args.pop('module')
        return module(**args)
