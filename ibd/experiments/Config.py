import yaml


class Config:
    def __init__(self, name, directory):
        self.name = name
        self.directory = directory

    def get_value(self, key):
        return self.config.get(key)

    def set_value(self, key, value):
        self.config[key] = value

    def save_config(self):
        with open(self.filepath, 'w') as file:
            yaml.dump(self.config, file)

    @staticmethod
    def load_config(name, directory):
        with open(f'{directory}/{name}.yaml', 'r') as file:
            config = yaml.safe_load(file)
        return config
