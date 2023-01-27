---
title: Tesla Dojo, the Race for AGI Dominance
author: gesit
date: 2023-01-26 14:10:00 +0800
categories: [blogging, ai]
tags: [ai, hardware] #has to be lowercase
render_with_liquid: false
---

A speculation on how compute latency might determine who will dominate upcoming AGI power struggle. 

## What Will AGI Look Like 
Jim Keller, the architect of many computing technology over the course of 40 years, think that the current neural network based AI can be viewed as a very large parallel system with sparse global memory. 


## How Will GPT like Model Run in Dojo 

[look at cerebras IO, booth has primary SRAM memory but vastly different IO bandwidth]

The first generation Dojo Exapod will have 1.1 TB of SRAM memory. Those memory are divided into chunk the size of 1.25 MB that reside inside a single core custom RISC-V like CPU. The CPU can churn 1 TOPS[footnote] bf16 with the ability to read 400 GB/s and write 270 GB/s those SRAM. The SRAM has 8 and 16 Byte of granularity.

Next, 354 of those cores are arranged into 19 X ? grid and integrated together to form a single die with a massive mesh network. The network communicates with the SRAM at 64 Byte per cycle in each read/write direction, so at 2 GHz a Dojo core can feed and produce 128 GB of data per second. The network mesh run bi directional trough every core and capable of 128 Byte pe cycle, so in theory 128*4*2 bilion = 1 TB of unique data passes trough every core. 
The D1 Chip measure only about 3 cm at its diagonal, so with respect to the speed of light at copper and 2GHz ferquency those communication can be done within one clock cycle or probably lose to it.






It is at this point I realize that these guys are building this computer to win AGI.



















































