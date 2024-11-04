import argparse

parser = argparse.ArgumentParser(prog='Test Program',description='This is a test',)
parser.add_argument('--input', required=True,help='This a test input')

args = parser.parse_args()

print(f'I: {args.input}')