# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 19:57:13 2021

@author: mehdi
"""
import matplotlib.pyplot as plt

figure, axes = plt.subplots()
#plt.suptitle(name_instn)
#for o in range(len(B)):
#    plt.scatter(B[o][0], B[o][1], c = "blue")
plt.scatter(3.0475483, 5.95591749, marker = "X", c = 'red')
plt.scatter(4.21722528, 7.25261494, marker = "X", c = 'orange')
plt.scatter(8.4424491, 0.53376056, marker = "X", c = 'orange')
plt.scatter(9.11330044, 1.67671293, marker = "X", c = 'red')
plt.plot((0,10), (11.28009217,-0.995000014), c = 'blue')

draw_circle0 = plt.Circle((3.0475483, 5.95591749), 1, fill=False)
draw_circle = plt.Circle((4.21722528, 7.25261494), 1, fill=False)
draw_circle1 = plt.Circle((8.4424491, 0.53376056), 1, fill=False)
draw_circle2 = plt.Circle((9.11330044, 1.67671293), 1, fill=False)
axes.set_aspect(1)
axes.add_artist(draw_circle)
axes.add_artist(draw_circle1)
axes.add_artist(draw_circle0)
axes.add_artist(draw_circle2)

plt.legend()