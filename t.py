def string_to_int(s):
    # 将字符串转换为字节
    bytes_data = s.encode('utf-8')
    # 将字节转换为整数
    return int.from_bytes(bytes_data, byteorder='big')

def encode_string(s):
    # 将左括号 '(' 替换为 '0'
    s = s.replace('(', '0')
    # 将左方括号 '[' 替换为 '1'
    s = s.replace(')', '1')
    s = s.replace('[', '2')
    s = s.replace(']', '3')
    return s

def decode_string(s):
    s = s.replace('0', '(')
    # 将左方括号 '[' 替换为 '1'
    s = s.replace('1', ')')
    s = s.replace('2', '[')
    s = s.replace('3', ']')
    return s

# 示例
original_str = "This is a test string with (brackets) and [square brackets]."
encoded_str = encode_string(original_str)
decoded_str = decode_string(encoded_str)

print(f"Original String: {original_str}")
print(f"Encoded String: {encoded_str}")
print(f"Decoded String: {decoded_str}")

# 示例
s = "[H]OC([H])([H])O|[H][O][C]([H])([H])[O]"
