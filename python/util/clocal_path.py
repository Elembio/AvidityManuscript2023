from pathlib import Path
from s3path import S3Path

# Abstraction of Cloud (S3) and local path
class ClocalPath:

    @classmethod
    def construct_path(cls, path):
        path = cls.as_string(path)
        if path is None:
            return None
        if path.startswith('s3://'):
            return S3Path(path[4:])
        elif path.startswith('s3:/'):
            return S3Path(path[3:])
        else:
            return Path(path)

    @staticmethod
    def as_string(path):
        if path is None:
            return None
        elif issubclass(type(path), S3Path):
            return path.as_uri()
        elif issubclass(type(path), Path):
            return path.as_posix()
        elif type(path) is str:
            return path
        else:
            raise RuntimeError(f'cannot get path string from object which is not S3Path, Path or string type')
