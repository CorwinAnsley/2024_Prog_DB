import re
TEST_STR = '$17£29€8¥55'

REGEX = '(\W)(\d+)'

# 0.82 dollars to the pound
# 0.87 euros to the pound
# 0.11 yuan to the pound

EXCHANGE_DICT = {'$':0.82,
                 '€':0.87,
                 '¥':0.11,
                 '£':1}


def test_regex_iter(regex, str):
    results = re.finditer(regex, str)
    for r in results:
        print(r)

def convert_prices(str, regex, exchange_dict):
    matches = re.finditer(regex, str)
    for match in matches:
        print(match.group(0))
        pound_value = exchange_dict[match.group(1)]*int(match.group(2))
        print(f'Value in pounds: {pound_value}')

#test_regex_iter(REGEX,TEST_STR)
convert_prices(TEST_STR,REGEX,EXCHANGE_DICT)