**Universidade Federal do Ceará (UFC)** - http://www.ufc.br

**Mestrado e Doutorado em Ciência da Computação (MDCC)** - http://mdcc.ufc.br

**Grupo de Pesquisa em Paralelismo, Grafos e Otimização (ParGO)** - http://www.lia.ufc.br/~pargo/

# Aims and Scope #

This project aims to create versions of programs from the NAS Parallel Benchmarks (NPB) for evaluating the performance of [HPE (Hash Programming Environment)](http://code.google.com/p/hash-programming-environment/), a platform of parallel components targeted at cluster computing platforms, developed at the Federal Universitity of Ceará (UFC) by the ParGO research group.

The project will be executed in the following stages:

  1. To port JAVA versions of the NPB programs to C#, since HPE is implemented on top of the Mono platform.
  1. To paralellize C# versions using MPI.NET. For that, the parallelization strategies adopted in MPI versions (C and Fortran) will be implemented.
  1. To refactor the MPI.NET versions in #-components and deploy them to an HPE Back-End for performance evaluation.

Besides to evaluate performance of HPE, it is also intended to make comparative studies about the performance of virtual execution platforms based on the JVM and CLI standards, in the context of applications with HPC requirements.

Also, it will be possible to evaluate performance of MPI.NET, the .NET/Mono binding for MPI (http://www.osl.iu.edu/research/mpi.net/) that has been used to implement message passing synchronizer in HPE.

# Achieved Results (Partial) #

  * Parallel versions of **SP**, **BT**, **LU** and **FT** in C#, using [MPI.NET](http://osl.iu.edu/research/mpi.net/) for enabling parallelism, have been built. The parallelism strategy adopted is the same one applied in their Fortran/MPI versions. The code structure is not object-oriented. It is a monolithic and procedural code written in C#, similar to the Fortran version, avoiding any performance overhead caused by an object-oriented structure;
    * **SP**: http://npb-for-hpe.googlecode.com/svn/trunk/NPB-CSharp/SP.MPI.HASH/
    * **BT**: http://npb-for-hpe.googlecode.com/svn/trunk/NPB-CSharp/BT.MPI.HASH/
    * **LU**: http://npb-for-hpe.googlecode.com/svn/trunk/NPB-CSharp/LU.MPI.HASH/
    * **FT**: http://npb-for-hpe.googlecode.com/svn/trunk/NPB-CSharp/FT.MPI.HASH/

  * Versions of **BT**, **SP**, **LU** and **FT** have been refactored onto components for the [HPE platform](http://hash-programming-environment.googlecode.com). The source code of the components (configurations and C# code of units) can be found in the repository at http://npb-for-hpe.googlecode.com/svn/branches/hpe_version/. More information about the HPE versions can be retrieved from the wiki page at http://code.google.com/p/npb-for-hpe/wiki/WikiCOPPversions.

  * [Click here](https://docs.google.com/spreadsheet/ccc?key=0AqD_KVtCMyYAdGhWOUoweEQ4T2JaZ1c4WU9ubVpfd0E&hl=en_US) to see a spreadsheet containing all the data collected in the performance evaluation to compare the component-based and the monolithic versions of **SP**, **BT**, **LU**, and **FT**.

  * Performance evaluation of JVM and CLI virtual machine implementations for common approaches to implement multidimensional arrays. This work has been submitted to SBAC-PAD 2012. The source codes used in the experiment are available at http://npb-for-hpe.googlecode.com/svn/branches/experimento/. Here a description of the folders that can be found there:
    * **v1** - versions AU (unidimensional arrays with inline index calculation arithmetic), in Java and C#;
    * **v1.2** - versions AU`*` (unidimensional arrays with index calculation arithmetic separated in methods), in Java;
    * **v2** - versions AJC and AJR (row-major and column-major jagged arrays), in Java and C#;
    * **v3** - versions ARR and ARC (row-major and column-major rectangular arrays), in C#;
    * **scripts** - the scripts used to compile and run the programs (must be adapted for your needs).


