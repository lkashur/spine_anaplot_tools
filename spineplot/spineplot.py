from argparse import ArgumentParser
from analysis import Analysis

def main(config, input):
    ana = Analysis(config, input)
    ana.run()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--config", help="Path to the configuration file")
    parser.add_argument("--input", help="Path to the input ROOT file")

    args = parser.parse_args()
    main(args.config, args.input)