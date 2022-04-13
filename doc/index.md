# CPLibrary

This Library is composents of many modules that interacts with each other:

- [Algebra module](algebra/index.md)
- [Linear Algebra module](linear_algebra/index.md)
- [Functional module](functional/index.md)
- [Machine Learning module]()
- [Number Theory module](nt/index.md)
- Polynomial module
- Topology module

```mermaid
flowchart BT
    Alg["Algebra"]    	
    LinAlg["Linear Algebra"]	
    Poly["Polynomial"]-->|Factoriser:FFT|Num
	Functional
    Num["Number Theory"]
    Poly & Num-->|egcd|Alg
   	LinAlg-->|Polynomial:Characteristic Polynomial|Poly
   	Topo["Topology"]-->LinAlg
   	Num-->|Ring Extension:Modular Sqrt|Poly
    DS["Data Structures"]-->|Order:Order Statistic Tree\nBinary Operation|Alg
    ML["Machine Learning"]-->Topo & LinAlg & Functional

```



