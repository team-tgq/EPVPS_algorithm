from datetime import datetime

a = datetime.now()

for i in range(100000000):
    b = 1

b = datetime.now()
c = b - a
print(c)
