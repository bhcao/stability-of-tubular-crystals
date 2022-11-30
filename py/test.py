import nano

m = 13
n = 13
repeat = 4
direction = 3
glide = 2
climb = 2
bn = 0
nano.init(m, n, repeat, direction, glide, climb, bn)

nano.dump("test", 0b111000111)

nano.delete()