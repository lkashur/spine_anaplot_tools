# imports
import numpy as np
import uproot
import argparse

def main(args):

    rf = uproot.open(args.input)
    
    pot_mc = rf['events/mc/POT'].to_numpy()[0][0]
    livetime_mc = rf['events/mc/Livetime'].to_numpy()[0][0]

    pot_onbeam = rf['events/onbeam/POT'].to_numpy()[0][0]
    livetime_onbeam = rf['events/onbeam/Livetime'].to_numpy()[0][0]

    pot_offbeam = rf['events/offbeam/POT'].to_numpy()[0][0]
    livetime_offbeam = rf['events/offbeam/Livetime'].to_numpy()[0][0]

    print(f'MC POT: {pot_mc}')
    print(f'MC Livetime: {livetime_mc}')
    print(f'On-beam POT: {pot_onbeam}')
    print(f'On-beam Livetime: {livetime_onbeam}')
    print(f'Off-beam POT: {pot_offbeam}')
    print(f'Off-beam Livetime: {livetime_offbeam}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    args = parser.parse_args()
    main(args)
