from argparse import ArgumentParser
from analysis import Analysis

def main(config, input):
    ana = Analysis(config, input)
    ana.override_exposure("offbeam", 534163 * (0.0309638 / 0.0240737) * (1.92082e19 / 159049), exposure_type='pot')
    ana.run()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--config", help="Path to the configuration file")
    parser.add_argument("--input", help="Path to the input ROOT file")

    args = parser.parse_args()
    main(args.config, args.input)