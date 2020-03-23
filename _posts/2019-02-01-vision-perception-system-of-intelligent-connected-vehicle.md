---
title: Vision Perception System of Intelligent Connected Vehicle
categories:
- Project
- Intelligent Connected Vehicle
---

From Setempber 1st, 2018 to January 31st, 2019, I participated in our college project of Intelligent Connected Vehicle (ICV). In this project, I took charge of the vision perception system of the ICV. I finished the whole perception system setup and calibration of multi-sensors. Our recording platform is an electric vehicle, which has been modified with actuators for pedals (acceleration and brake) and steering wheel. We use the following sensors:

- One Laser scanner: [RS-LIDAR-32](http://www.robosense.cn/rslidar/rs-lidar-32)
- Two Color cameras, 3.1 Megapixels: [MER-310-12UC](http://www.daheng-image.com/products/ProductDetails.aspx?current=5&productid=2649)

The laser scanner spins at 10 frames per second, capturing approximately 600k points per second. The vertical resolution of the laser scanner is 32. The cameras are mounted approximately level with the ground plane. The cameras are triggered at 10 frames per second by the laser scanner (when facing forward) with shutter time adjusted dynamically (maximum shutter time: 2ms).

{% include figure.html image="/assets/imgs/sdc1.jpg" %}

{% include figure.html image="/assets/imgs/sdc2.jpg" %}

{% include figure.html image="/assets/imgs/sdc3.jpg" %}