string = "GCNAT"

if 'N' in string:
    string_list = list(string)
    string_list.remove('N')
    string_list.append('A')

print(''.join(string_list))