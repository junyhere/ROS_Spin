"""Legacy compatibility boundary; not part of the active chemical model."""
from legacy_entrypoint import retired_entrypoint


def main() -> None:
    retired_entrypoint("ROS_Util.py")


if __name__ == "__main__":
    main()
