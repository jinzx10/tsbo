import re


pattern  = re.compile(r'([^ ]+)([ \t]+)([^#]*[^ #]|)(#.*)?')

line = 'key    value   # comment'
keyval = re.search(pattern, line).group(1,3)
print(keyval)

line = 'key    # comment'
keyval = re.search(pattern, line).group(1,3)
print(keyval)

line = 'key    '
keyval = re.search(pattern, line).group(1,3)
print(keyval)

line = 'key    3 3 3 # comment'
keyval = re.search(pattern, line).group(1,3)
print(keyval)

