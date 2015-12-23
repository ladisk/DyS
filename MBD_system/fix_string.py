__author__ = 'lskrinjar'

def fix_string(_str):
    """

    :param _str:
    :return:
    """
    _str = _str.title()

    if "Mc" in _str:
        #   find character to be replaced
        char_to_be_replaced = _str[_str.index('Mc')+2]
        #   replace lower case character with upper case character
        _str = _str.replace(char_to_be_replaced, char_to_be_replaced.upper())

    elif "Mac" in _str:
        #   find character to be replaced
        char_to_be_replaced =  _str[_str.index('Mac')+3]
        #   replace lower case character with upper case character
        _str = _str.replace(char_to_be_replaced, char_to_be_replaced.upper())

    if "Et" in _str:
        #   find character to be replaced
        char_to_be_replaced =  _str[_str.index('Et')]
        #   replace lower case character with upper case character
        _str = _str.replace(char_to_be_replaced, char_to_be_replaced.lower())

    if "Al" in _str:
        #   find character to be replaced
        char_to_be_replaced =  _str[_str.index('Al')]
        #   replace lower case character with upper case character
        _str = _str.replace(char_to_be_replaced, char_to_be_replaced.lower())

    return _str

if __name__ == "__main__":
    # print fix_string("Herbert-Mcwhannell")
    print fix_string("Flores Et Al")