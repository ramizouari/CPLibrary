# CPLibrary

This Library is composents of many modules that interacts with each other:

- [Algebra module](linear_algebra/index.md)
- Linear Algebra module
- Functional module
- Machine Learning module
- Number Theory module
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



```mermaid
flowchart LR
    subgraph Alg["Algebra"]
    	subgraph Abstract
    	subgraph Integral Domain
    	end
    	end
    	subgraph Order
    	O1["Order Closure"]
    	O2["Order Operations"]
    	end
    	subgraph Binary Operation
    	end
    end
    subgraph LinAlg["Linear Algebra"]
        subgraph Vec["Vector Space"]
        end
        subgraph Mat["Matrix"]
        	CharPoly["Characteristic Poly"] -->Polynomial
        end
    end
    subgraph Functional
    end
    subgraph Number Theory
    	subgraph Factorisation
    	Factoriser
    	end
    	subgraph Modular Arithmetic
    	UnityRoot["Primitve Root"]-->Factoriser
        SQRT["Square Root"] --> QuadraticExtension
    	end
    end
    subgraph Data Structures
    OST["Order Statistics Tree"] -->O1
    end
    subgraph Poly["Polynomial"]
    	subgraph Polynomial
    	end
    	subgraph Ring Extension
    	end
    end

```