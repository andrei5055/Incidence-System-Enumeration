<TeXmacs|2.1>

<style|<tuple|article|number-long-article>>

<\body>
  <doc-data|<doc-title|Enumeration of Combined
  BIBDs>|<doc-subtitle|>|<doc-author|<author-data|<author-name|Andrei
  Ivanov>|<\author-email>
    <em|andrey.ivanov@gmail.com>
  </author-email>>>>

  <abstract-data|<\abstract>
    A BIBD (Balanced Incomplete Block Design)
    <math|D=<around*|(|V,\<bbb-B\>|)> is called<embold|<text|<em|combined>>>>,
    if <math|for some<hspace|2> <with|math-display|true|<line-break>n\<geqslant\><no-break><with|math-condensed|false|>2>
    the set of blocks \<bbb-B\> is <text|an a>disjunctive association of
    <around*|{|\<bbb-B\><rsub|i>|}> such that
    every<line-break>D<rsub|i>=<around*|(|V,<no-break>\<bbb-B\><rsub|i>|)> is
    <text|a> BIBD.<text| In this paper we will describe the algorithm of
    enumeration of combined BIBDs and some of the results, obtained with
    it.><rsub|\<nosymbol\>>>

    \;
  </abstract>>

  <section|Introduction>

  Let <math|V> be a set of the elements and \<bbb-B\> is a collection of some
  (not necessary different) subsets of <math|V.> We will denote by
  <math|v=<around|\||V|\|>> and <math|b=<around|\||\<bbb-B\>|\|>> the
  cardinalities of these sets. Usually, the elements
  <math|B<rsub|j>\<in\>\<bbb-B\>,1\<leqslant\>j\<leqslant\>b> are called the
  blocks and the pair <math|D>=(<math|V>, <math|\<bbb-B\>>) is called a
  (combinatorial) <em|<strong|block design>> or an
  <em|<strong|<em|<em|i<em|<em|ncidence system>>>><em|>>>.\ 

  \;

  <\with|par-first|0tab>
    For a block design <math|D>=(<math|V>, <math|\<bbb-B\>>) the number of
    blocks that have the same set of elements as the block
    <math|B<rsub|\<nosymbol\>>\<in\>\<bbb-B\> will be called> <strong|<em|the
    multiplicity>> of the that block:
  </with>

  <\equation*>
    \<mu\><around*|(|B<rsub|\<nosymbol\>>|)>=<around*|\||<around*|{||\<nobracket\>>*B<rsub|i>\<in\>\<bbb-B\>:B<rsub|i>=B<rsub|\<nosymbol\>>|}><around*|\||.|\<nobracket\>>
  </equation*>

  <\with|par-first|0tab>
    Two block designs <math|D<rsub|i>>=(<math|V<rsub|i>,\<bbb-B\><rsub|i>>),
    <math|i = 1,2> <math|>are said to be <strong|<em|isomorphic>> if
    <math|<around*|\||V<rsub|1>|\|>=<around*|\||V<rsub|2>|\|>> and there
    exists a bijection function <math|\<frak-I\>:V<rsub|1>\<rightarrow\>V<rsub|2>>
    such that for <math|any element x\<in\>V<rsub|1> and any block
    B<rsub|1>\<in\>\<bbb-B\><rsub|1>> of size
    <math|k<rsub|1>=<around*|\||B<rsub|1>|\|>> which contains <math|x>, the
    element <math|\<frak-I\><around*|(|x|)> belongs exactly to
    \<mu\><around*|(|B<rsub|1>|)> blocks B<rsub|2>\<in\>\<bbb-B\><rsub|2>> of
    the same size <math|k<rsub|1>>.<space|1em>

    <\equation*>
      <around*|{||\<nobracket\>><around*|{|\<frak-I\><around*|(|x|)>:x\<in\>B<rsub|1><text|}>:B<rsub|1>\<in\>\<bbb-B\><rsub|1><text|}><rsub|\<nosymbol\>>
      \<mu\><around*|(|B<rsub|1><rsub|\<nosymbol\>>|)>=\<mu\><around*|(|B<rsub|2>|)>|}>=\<bbb-B\><rsub|2.>
    </equation*>

    In other words, if we rename every element <em|x> by
    <math|\<frak-I\>>(<em|x>), then the collection of blocks
    <math|\<bbb-B\><rsub|1>> is transformed into <math|\<bbb-B\><rsub|2>>.The
    bijection <math|\<frak-I\>> is called an <strong|<em|isomorphism>>.

    <strong|<strong|>>
  </with>

  <\with|par-first|0tab>
    Let \ <math|D>=(<math|V>, <math|\<bbb-B\>>) is a block design with
    <math|v=<around*|\||V|\|> elements<infix-and>b=<around*|\||\<bbb-B\>|\|><text|>>
    blocks. The <strong|<em|<strong|incidence matrix>><em|>> of this block
    design is the<nbsp><math|v \<times\>b> matrix <math|M(D)> whose entries
    are defined as
  </with>

  <\equation>
    m<rsub|i j>=<choice|<tformat|<table|<row|<cell|1, if
    x<rsub|i>\<in\>B<rsub|j><space|2em>>>|<row|<cell|0, otherwise.>>>>>
  </equation>

  <\with|par-first|0tab>
    \;

    It is easy to prove that two block designs
    <math|D<rsub|i>>=(<math|V<rsub|i>,\<bbb-B\><rsub|i>>), <math|i = 1,2>,
    <math|>are isomorphic if and only if there are such permutation matrices
    <math|P<infix-and>Q, acting,respectively, on the rows of <infix-and>
    columns of \ incidence matrices M<around*|(|D<rsub|i>|)> that the
    following holds: >

    <\equation>
      P\<cdot\>M<around*|(|D<rsub|1>|)>=M<around*|(|D<rsub|2>|)> \<cdot\>Q.
    </equation>
  </with>

  <\with|par-first|0tab>
    Let <math|><math|D<rsub|<math|i>>=<around*|(|V,\<bbb-B\><rsub|<math|i>>|)>>,
    <math|i=1,2,><math|be two >block designs defined on the same set of
    elements<math| V> and <math|M<rsub|i>>=<math|M<rsub|i><around*|(|D<rsub|i>|)>
    are corresponding \ incidence matrices>. By horizontally joining
    (concatenating) these matrices we can construct
    <math|<around*|\||V|\|>\<times\><around*|(|<around*|\||\<bbb-B\><rsub|<math|1>>|\|>+<around*|\||\<bbb-B\><rsub|<math|2>>|\|>|)>>
    matrix <math|M<around*|(|D|)>=<around*|[|M<rsub|1>,M<rsub|2>|]>>,
    <math|<math|<math|which will be<text| an incidence >matrix > >>of some
    block design <math|D<rsub|<math|>>=<around*|(|V,\<bbb-B\><rsub|<math|1>>+\<bbb-B\><rsub|2>|)>>.\ 

    \;

    The block design (resp., incidence matrix) constructed by using such
    horizontal concatenation of the insidence matrices of the smaller block
    designs will be called <em|<em|<em|<strong|combined<em|>>><strong|>>>
    block design (resp., incidence matrix).
  </with>

  \;

  <\with|par-first|0tab>
    The block design <math|D>=(<math|V>, <math|\<bbb-B\>>) is called a<em|
    <strong|>><strong|<em|<strong|>balanced incomplete block design>> (BIBD),
    if
  </with>

  <\vgroup>
    (a) each block <math|B<rsub|j>\<in\>\<bbb-B\>> contains the same number
    of elements <em|k>;

    (b) each element <math|v<rsub|i>\<in\>V> belongs to the same number of
    blocks <math|r>;

    (c) any two of distinct elements <math|v<rsub|i>,v<rsub|j>\<in\>V> appear
    together in the same number of blocks <math|\<lambda\>>.
  </vgroup>

  \;

  <\with|par-first|0tab>
    These five parameters (<math|v,b,r,k,\<lambda\>>) satisfy following two
    equations:
  </with>

  <\equation*>
    v r=<rigid|>b k,
  </equation*>

  <\equation*>
    r <around*|(|k-1|)>=\<lambda\> <around*|(|v-1|)>,
  </equation*>

  which allow to express <math|r<infix-and>b> as

  <\equation>
    r=\<lambda\> <frac|v-1|k-1>,<text|<space|3em>>b=\<lambda\> <frac|v
    <around*|(|v-1|)>|k<around*|(|k-1|)>>
  </equation>

  Much more information regarding the BIBDs and related structures could be
  found in [1].

  \;

  <\with|par-first|0tab>
    Suppose, we do have two (not nessesary different or non-isomorfic)
    BIBD's: D<rsub|<math|i>>=<around*|(|V,\<bbb-B\><rsub|<math|i>>|)>,
    <math|i = 1,2> <math|,>with parameters
    <math|<around*|(|v,k,\<lambda\><rsub|i>|)>>. It means that not only all
    these designs are defined on the same set of elements <math|V,but also
    that <around*|\||X|\|>=<around*|\||Y|\|>=k,for any two block
    X\<in\>\<bbb-B\><rsub|1>,Y\<in\>\<bbb-B\><rsub|2>>.
  </with>

  \ 

  <\with|par-first|0tab>
    Let's choose some permutation matrices
    <math|P<rsub|1><infix-and>P<rsub|2>> acting respectively on the rows of
    incidence matrices <math|M<rsub|i>=M<rsub|><around*|(|D<rsub|i>|)>. >It
    is easy to check that the combined incidence matrix
    <math|M=<around*|[|P<rsub|1>\<cdot\>M<rsub|1>,P<rsub|2>\<cdot\>M<rsub|2>|]>>
    satisfies the following conditions:

    <\enumerate-alpha>
      <item>each column of <math|M> contains the exactly <em|k> 1's;

      <item>each row of <math|M> contains exactly
      <math|><math|r=r<rsub|1>+r<rsub|2>=<around*|(|\<lambda\><rsub|1>+\<lambda\><rsub|2>|)>
      <around*|(|v-1|)>/<around*|(|k-1|)>> \ 1's

      <item>two distinct rows of <math|M> contain both 1's in exactly
      <math|\<lambda\>=\<lambda\><rsub|1>+\<lambda\><rsub|2><space|1em>columns.>
    </enumerate-alpha>

    Therefore, the combined incidence matrix
    <math|M=<around*|[|P<rsub|1>\<cdot\>M<rsub|1>,P<rsub|2>\<cdot\>M<rsub|2>|]>>
    corresponds to some combined BIBD (CBIBD).

    \;

    The concatenation of incidence matrices of BIBD's having the same
    parameters (<math|v,k>) is probably the easiest way to construct bigger
    BIBD's from smaller ones. Using\ 

    <\enumerate-alpha>
      <item>different pairs <math|<around*|(|<wide|P|~><rsub|1>,
      <wide|P|~><rsub|2 >|)>, >generally speaking, we can construct pairwise
      non-isomorphic BIBDs with the same parameters;

      <item>two, three or more \Pinitial\Q BIBDs as units, we can construct
      BIBDs that are composed of two, three or more BIBD's of smaller sizes.
    </enumerate-alpha>
  </with>

  <paragraph|>

  <\with|par-first|0tab>
    A typical problem in combinatorics is the construction of a complete set
    of pairwise non-isomorphic combinatorial objects with given properties.
    This is commonly referred to as the problem of enumerating such objects.
    Our goal: to construct an algorithm that would allow one to enumerate
    CBIBDs for <math|n\<gtr\><no-break>1> and given set of parameters
    <math|(v, k, <around*|{|\<lambda\><rsub|1>,\<ldots\>,\<lambda\><rsub|n>|}>>),
    where <math|<around*|(|v, k, \<lambda\><rsub|i>|)> are the parameters of
    BIBD \ for 1\<leqslant\><no-break>i\<leqslant\><no-break>n.>
  </with>

  <\section>
    Representation of the Combined BIBDs
  </section>

  We will follow the methodology described in [2] and [5], where it is
  recommended to build a complete list of so-called <strong|<em|canonical>
  >objects for solving a particular combinatorial enumeration problem. Here
  is a brief description of this technique.

  \;

  <\with|par-first|0tab>
    Suppose we need to enumerate the combinatorial objects with some set of
    specified properties. Let \<cal-A\> be a set of all such objects and
    <math|G=<no-break>G<around*|(|\<cal-A\>|)>> be a group of transformations
    applied to objects from \<cal-A\> such that
    <math|g<around*|(|\<alpha\>|)> \<in\>\<cal-A\> >for
    <math|\<forall\>\<alpha\>\<in\><no-break>\<cal-A\><infix-and>><math|\<forall\>g
    \<in\>G.<text| Usually such a transformation group G is called
    <em|<strong|the isomorphism group>>>. >
  </with>

  \;

  <\with|par-first|0tab>
    In most cases, the isomorphism group <math|G> is somehow related to the
    renumbering of some components of objects from \<cal-A\>, <math|but it
    can also include transformations of some other types.>
  </with>

  \;

  <\with|par-first|0tab>
    Examples of transformations of the first type are the renumbering of
    vertices of ordinary graphs or elements and blocks of incidence systems,
    etc. Other examples of transformations that can also appear in <math|G>
    when enumerating combinatorial objects of the corresponding classes are:
    replacing edges with non-edges (and vice versa) in regular graphs of
    degree <math|k> with <math|v = 2k + 1> vertices, or replacing incidences
    with non-incidences (and vice versa) in BIBDs with <math|v = 2k>
    elements, or even the transposition of the incidence matrices of
    symmetric BIBDs.
  </with>

  \;

  <\with|par-first|0tab>
    Suppose that each object <math|\<alpha\>>\<in\>\<cal-A\> has some matrix
    representation <math|M<around*|(|\<alpha\>|)>>. It could be an adjacency
    matrix of a graph, an incidence matrix of a block design, an
    <math|n\<times\>n> matrix representing <text|a> Latin square etc. For
    matrices <math|M<around*|(|\<alpha\>|)>,<text| >M<around*|(|\<beta\>|)>>
    \ we can determine the binary precedence relation, for example, by
    lexicographic comparison of their rows. We will say that <math|the matrix
    M<around*|(|\<alpha\>|)> is <em|<text|<strong|greater<em|>><strong|>>>
    than the matrix M<around*|(|\<beta\>|)>> if for some
    <math|n<math|\<geqslant\>0>> their first<math| <em|n>> rows are
    identical, but (n+1)-th row of <math|M<around*|(|\<alpha\>|)> is
    lexicographically greater than >(n+1)-th row of
    <math|M<around*|(|\<beta\>|)>. This fact will be expressed as:>

    <\equation*>
      M<around*|(|\<alpha\>|)>\<succ\>M<around*|(|\<beta\>|)>
    </equation*>
  </with>

  <\with|par-first|0tab>
    <em|<em|<strong|Definition 2.1: > <strong|>><strong|>>A combinatorial
    object <math|\<alpha\>>\<in\>\<cal-A\> (and its matrix representation
    <math|M<around*|(|\<alpha\>|)>>) will be called <strong|<em|canonical>>
    with respect to the isomorphism group <em|G>, if for <math|\<forall\>g
    \<in\>G>:
  </with>

  <\equation*>
    M<around*|(|\<alpha\>|)>\<succcurlyeq\>M<around*|(|g<around*|(|\<alpha\>|)>|)>.
  </equation*>

  By definition, the canonicity or non-canonicity of the combinatorial oject
  <math|\<alpha\>\<in\> >\<cal-A\> could be detrnined by its matrix
  representation <math|M<around*|(|\<alpha\>|)><infix-and>the >group
  <math|G>(\<cal-A\>). In many cases, this fact greatly simplifies the
  construction of all existing pairwise non-isomorphic combinatorial objects
  of the class \<cal-A\>, since there is no need to maintain lists of already
  constructed objects and check for isomorphism for each newly constructed
  object.

  \;

  <\with|par-first|0tab>
    <em|<em|<strong|Definition 2.2: > <strong|>><strong|>>A balanced
    incomplete block design <math|D = (V,\<bbb-B\>)> will be called
    <em|<strong|combined of rank n> <strong|>>if the set of block
    <math|<em|>\<bbb-B\>> can be divided into <math|n\<geqslant\>2> disjoint
    subsets <math|<around*|(|components|)> <around*|{|\<bbb-B\><rsub|i>|}>>,
    in such a way that <math| D<rsub|i>=(V,<no-break><text|<no-break>>\<bbb-B\><rsub|i>)
    will be <text| a> balanced incomplete block design for each
    1\<leqslant\>i\<leqslant\>n>.
  </with>

  \;

  <\with|par-first|0tab>
    Since <math|D<rsub|i>=(V,<no-break><text|<no-break>>\<bbb-B\><rsub|i>) is
    <text|a> BIBD with parameters <around*|(|v,k,\<lambda\><lsub|i>|)> >for
    each <math|0\<leqslant\>i\<leqslant\>n>, it's very logical to use
    <no-indent>

    <\equation>
      (<math|v, k,<around*|{|\<lambda\><lsub|1>,\<lambda\><lsub|2>,\<ldots\>,\<lambda\><lsub|n>|}>>)
    </equation>

    \ as the parameters of combined BIBD <math|D = (V,\<bbb-B\>)>.\ 
  </with>

  \;

  <\with|par-first|0tab>
    We will represent a combined BIBD <em|D> by the the pair
    <math|<around*|(|\<frak-m\>,M<around*|(|D|)>|)>><math|, \ where
    \<frak-m\>>=(<math|\<frak-m\><rsub|j>>) is a vector of length
    \|<math|\<bbb-B\><around*|\||,with integer coordinates
    |\<nobracket\>>>indicating the occurrence of corresponding block in
    <math|\<bbb-B\><rsub|i>,M<around*|(|D|)> is the incidence matrix of D
    .<text|>>
  </with>

  \;

  <\with|par-first|0tab>
    In what follows, we will represent this pair
    <math|<around*|(|m,M<around*|(|D|)>|)>> in one matrix:
  </with>

  <\equation>
    \<bbb-M\><around*|(|D|)>=<around*|[|<with|font-base-size|12|<below|<with|font-base-size|14|<above||<with|font-base-size|14|<with|font-base-size|14|m<rsub|>>>>>|<with|font-base-size|12|M<around*|(|D|)>>>>|]><space|1em>where
    m=<around*|(|m<rsub|j>|)><infix-and>m<rsub|j>= i, if B<rsub|j>
    \<in\>\<bbb-B\><rsub|i> for 1\<leqslant\>j\<leqslant\><around*|\||\<bbb-B\>|\|>
  </equation>

  \;

  <\with|par-first|0tab>
    with \|V\|+1 rows and \|\<bbb-B\>\| columns, and we will call (2.1)
    <em|<strong|the><strong|>> <em|<strong|matrix representation>> of <em|D>.
  </with>

  \;

  <\with|par-first|0tab>
    To formulate the CBIBD enumeration problem precisely, we also need
  </with>

  a) determine the group of isomorphisms acting on these objects;

  b) clarify some parameters of the enumerated CBIBDs.

  \;

  <\with|par-first|0tab>
    Let us first focus on (a) and consider the direct product
    <math|\<frak-G\>=><verbatim|<math|\<cal-V\> \<times\>\<cal-B\> >>of
    symmetric groups <verbatim|<math|\<cal-V\> >>and <math|\<cal-B\>>, acting
    on the rows and columns of matrix <math|M<around*|(|D|)>>, respectively.
    In accordance with the Definition 2.1, the matrix
    <math|\<bbb-M\><around*|(|D|)> will be called canonical with respect to
    the >isomorphism group <math|\<frak-G\>>, if\ 
  </with>

  <\equation*>
    \<bbb-M\><around*|(|D|)> \<succcurlyeq\>\<frak-g\><around*|(||\<nobracket\>>\<bbb-M\><around*|(|D|)><around*|)|,for
    |\<nobracket\>>any \<frak-g\>\<in\>\<frak-G\>.
  </equation*>

  <\proposition>
    If <math|\<bbb-M\><around*|(|D|)> <text|defined by
    <em|<around*|(|2.1|)>>> \ <em|>>is a canonical matrix of CBIBD of rank n
    then the coordinates of the vector <math|m=<around*|(|m<rsub|j>|)>> are\ 

    <\equation>
      <around*|(|m<rsub|j>|)>=<around*|(|n,\<ldots\>,n,n-1,\<ldots\>,n-1,n-2,\<ldots\>n-2,\<ldots\>,2,\<ldots\>,2,1,\<ldots\>1|)>
    </equation>

    which means that all blocks of<math|<normal-size|<very-large|>> >BIBD
    <math|D<rsub|i>=(V,<no-break><text|<no-break>>\<bbb-B\><rsub|i>)> are
    represented in <math|M<around*|(|D|)> <text|by the columns with
    consecutive indices<em|>> > <math|<text|for any
    <math|1\<leqslant\>i\<leqslant\>n. >>>
  </proposition>

  <\with|par-mode|left|par-first|0tab>
    <with|par-first|0tab|<em|Proof:> Let's consider the set of indices<math|
    I<rsub|n>=>{<math|i<rsub|1>,i<rsub|2>,<text-dots>,><em|i<math|<rsub|<around*|\||\<bbb-B\><rsub|n>|\|>>>>}
    of the columns of <math|\<bbb-M\><around*|(|D|)> >which correspond to all
    blocks from <math|<no-break><text|<no-break>>\<bbb-B\><rsub|n >>. By
    definition of<math| vector \ m=<around*|(|m<rsub|j>|)>>, the equality
    <math|m<rsub|l>=n is true for \<forall\>l<text|>\<in\><text|<no-indent>><no-break>I<rsub|n><text|<no-break>.>>>
  </with>

  \;

  <\with|par-first|0tab>
    Since <math|n is <text|a<em|>> maximal value among all >coordinates
    <math|<around*|(|m<rsub|j>|)>>, <math|1\<leqslant\>j\<leqslant\><around*|\||\<bbb-B\>|\|><infix-and>the
    >matrix <math|\<bbb-M\><around*|(|D|)> is canonical,<text|any index
    <em|l><em|<math|\<in\>>I><math|<rsub|n>>>> must precede any index
    <math|x> with a smaller value <math|m<rsub|x>\<less\>n.> Otherwise, by
    rearranging the columns<em| l> and <em|x> of the matrix \<bbb-M\>(<em|D>)
    we could get the matrix <math|<wide|\<bbb-M\>|~>>(<em|D>) with the first
    row <math|<wide|m |~>>\<succ\><no-break><em|m>, which contradicts our
    assumption about the canonicity of \<bbb-M\>(<em|D>). \ 

    \;

    Using the same arguments for \ <math|<no-break><text|<no-break>>\<bbb-B\><rsub|n-1>,><math|<no-break><text|<no-break>>\<bbb-B\><rsub|n-2>,\<ldots\>,<no-break><text|<no-break>>\<bbb-B\><rsub|1>,
    we >get (2.3).

    \;
  </with>

  <\with|par-first|0tab>
    Denote by <math|b<lsub|i>>=\|<math|<no-break><text|<no-break>>\<bbb-B\><rsub|i>>\|,
    the number of blocks of <math|i>-th component
    <math|D<rsub|i>=(V,<no-break><text|<no-break>>\<bbb-B\><rsub|i>),1\<leqslant\>i\<leqslant\>n
    of CBIBD.>

    \;

    <paragraph|Proposition 2.2.><em|If> <math|\<bbb-M\><around*|(|D|)>
    <text|<em|defined by> <em|<em|<around*|(|2.1|)>>>> \ <em|>><em|is a
    canonical matrix of CBIBD of rank > <math|n,>
    <em|the><em|n<space|1em>><math|b<math|<lsub|s
    >><text|<no-indent*>>\<geqslant\><text|<no-indent*>>b<lsub|t><lsup|>>
    <em|for any> <math|1\<leqslant\>s\<leqslant\><text|<math|t\<leqslant\>n.>>>

    \;

    <\with|par-mode|justify|par-first|0tab>
      <\with|par-mode|justify|par-sep|0fn|par-hyphen|normal>
        <text|<text|<with|math-display|true|<\with|math-condensed|false|par-sep|0.2fn|par-flexibility|1000|par-spacing|plain|par-hyphen|professional>
          <text|<text|<text|<em|Proof:> By definition,
          <math|b<lsub|i>>=\|<math|<no-break><text|<no-break>>\<bbb-B\><rsub|i>>,
          <math|1\<leqslant\>i\<leqslant\>n >is equal to the number of
          occurrences of <math|i> in the vector
          <math|m<text|=><around*|(|m<rsub|j>|)>. According to
          <rsub|>Proposition 2.1,<text|when <math|\<bbb-M\><around*|(|D|)> is
          canonical,>all coordinates >of this vector satisfy
          <around*|(|2.2|)>,which can be rewritten as<text|><text|> >>>>
        </with>>>>
      </with>

      <\equation*>
        <around*|(|m<rsub|j>|)>=<around*|(|<around*|[|b<lsub|1>**\<ast\>n|]>,\<ldots\>,<around*|[|b<lsub|p>\<ast\><around*|(|n-p+1|)>*,\<ldots\>,<around*|[|b<lsub|n>\<ast\>1|]>|)><text|>,|\<nobracket\>>
      </equation*>

      <math|where <around*|[|b<lsub|a>\<ast\><around*|(|n-a+1|)>|]> is
      <text|a> vector with all b<lsub|a> coordinates equal to
      <around*|(|n-a+1|)>.>

      \;

      <math|If the statement of >Proposition 2.2 is false, then for some
      <math|s> and <math| t>, \ <math|1\<leqslant\>s\<leqslant\><text|<math|t\<leqslant\>n,>>>
      <math|b<math|<lsub|s >><text|<no-indent*>>\<less\><text|<no-indent*>>b<lsub|t>>
      we do have

      <\equation*>
        <around*|(|m<rsub|j>|)>=<around*|(|<around*|[|b<lsub|1>**\<ast\>n|]>,\<ldots\>,<around*|[|b<lsub|s>\<ast\><around*|(|n-s+1|)>|]>*,\<ldots\>,<around*|[|b<lsub|t>\<ast\><around*|(|n-t+1|)>|]>,\<ldots\>.,<around*|[|b<lsub|n>\<ast\>1|]>|)><text|>.
      </equation*>

      <\padded-left-aligned>
        Keeping the indices of the components <math|{\<bbb-B\><lsub|i>},
        1\<leqslant\>i\<leqslant\>n, i\<neq\>s,t> the same and changing the
        indices of the components \<bbb-B\><math|<lsub|s>> and
        \<bbb-B\><math|<lsub|t>> by <math|t> and<math| s>, respectively, we
        will see that the vector
      </padded-left-aligned>

      <\equation*>
        <around*|(|<wide|m<rsub|j>|~>|)>=<around*|(|<around*|[|b<lsub|1>*\<ast\>n|]>,\<ldots\>,<around*|[|b<lsub|t>\<ast\><around*|(|n-s+1|)>|]>*,\<ldots\>,<around*|[|b<lsub|s>\<ast\><around*|(|n-t+1|)>|]>,\<ldots\>.,<around*|[|b<lsub|n>\<ast\>1|]>|)><text|>
      </equation*>

      is greater then <math|<around*|(|m<rsub|j>|)>>. For both of these
      vectors the first <math|<around*|(|n-s|)> blocks of their
      coordinates,equal to <text|><text|<no-indent>
      >n,<around*|(|n-1|)>,\<ldots\>,<around*|(|s+1|)>,respectively,>will be
      the same and the block <math|of the vector
      <around*|(|<wide|m<rsub|j>|~>|)>> with the coordinates equal to <math|s
      \ will be longer,since > <math|b<math|<lsub|t
      >><text|<no-indent*>>\<gtr\><text|<no-indent*>>b<lsub|s>>. Thus, we
      have obtained the relation <math|<around*|(|<wide|m<rsub|j>|~>|)>
      \<succ\>><math|<around*|(|m<rsub|j>|)>> which contradicts our
      assumption about the canonicity of the matrix
      \ <math|\<bbb-M\><around*|(|D|)>.>

      \;
    </with>

    A direct consequence of Proposition 2.2 is the following

    \ 

    <strong|Proposition 2.3.><em| For the canonical CBIBD of rank n with
    parameters> (<math|v, k,<around*|(|\<lambda\><lsub|1>,\<lambda\><lsub|2>,\<ldots\>\<lambda\><lsub|n>|)>>):

    <\equation*>
      \<lambda\><lsub|i>\<geqslant\>\<lambda\><lsub|j>
      <space|1em><em|<text|<text|for<em|>><em|>>><space|1em>\<forall\><around*|(|i,j|)>:1\<leqslant\>i\<leqslant\>j\<leqslant\>n.
    </equation*>

    <em|Proof:> According to (1.3)

    <\equation*>
      <text|>b<rsub|i>=\<lambda\><rsub|i> <frac|v
      <rsub|i><around*|(|v<rsub|i>-1|)>|k<rsub|i><around*|(|k<rsub|i>-1|)>>=<text|>\<lambda\><rsub|i>
      <frac|v <around*|(|v-1|)>|k<around*|(|k-1|)>><space|1em>for
      1\<leqslant\>i\<leqslant\>n.
    </equation*>

    As we could see, all <math|\<lambda\><lsub|i> are proportional to
    ><math|<text|>b<rsub|i>> with the same proportionality factor:
    <math|k<around*|(|k-1|)>/<around*|(|v<around*|(|v-<no-break>1|)>|)>.<text|
    But according to Proposition 2.2 <math|b<rsub|i>\<geqslant\>b<rsub|j >
    >for <math|\<forall\><around*|(|i,j|)>:1\<leqslant\>i\<leqslant\>j\<leqslant\>n><math|.
    Therefore,\<lambda\><lsub|i>\<geqslant\>\<lambda\><lsub|j>>.>>

    <\section>
      Isomorphisms of Combined BIBDs
    </section>

    If <em|D >is a canonical CBIBD with parameters \ (<math|v,
    k,<around*|{|\<lambda\><lsub|1>,\<lambda\><lsub|2>,\<ldots\>\<lambda\><lsub|n>|}>>),
    then, \ according to Propositions 2.1 and 2.2, its matrix representation
    must have the following form:

    <\eqnarray*>
      <tformat|<cwith|1|1|3|3|cell-halign|r>|<table|<row|<cell|\<bbb-M\><around*|(|D|)>=<around*|[|<with|font-base-size|12|<below|<with|font-base-size|14|<above||<with|font-base-size|14|<with|font-base-size|16|n\<ldots\>n<space|2em><space|1em><around*|(|n-1|)>\<ldots\><around*|(|n-1|)><space|1em>
      \<ldots\><space|2em>1. \ . . 1><rsub|>>>>|<with|font-base-size|18|<em|B<rsub|1><space|2em>
      ><space|3em><em|B<rsub|2><space|5em>\<ldots\><space|1em>><space|1em><em|B<rsub|<em|<em|n>>>
      >>>>|]>>|<cell|>|<cell|<text|> <around*|(|3.1|)>>>>>
    </eqnarray*>

    \;

    where <em|B<rsub|i><em|>> is a balanced incomplete block design with
    parameters \ (<math|v, k,\<lambda\><lsub|i>>) for
    <math|1\<leqslant\>i\<leqslant\>n.> Let's describe the group of
    transformations <math|G<around*|(|D|)> >which will be used for canonicity
    check of <em|D>.\ 

    \;

    First of all, as in the case of BIBD, it is natural to include in this
    group the symmetric group <em|<math|S><rsub|v>>=<math|S><rsub|<em|v>>(<em|V>),
    acting on the set of elements <em|V>.

    \;

    Second, <itemize|of course we should consider> the direct product of the
    groups: \ 

    <\equation*>
      <em|<text|S><around*|(|\<bbb-B\>|)>=>S<lsub|b<lsub|1>>\<times\>S<lsub|b<lsub|2>>\<times\>\<ldots\>\<times\>S<lsub|b<lsub|n>>,
    </equation*>

    where <math|S<lsub|b<lsub|i>>> = <math|S<lsub|b<lsub|i>>>(<math|B<lsub|i>>)
    <math|1\<leqslant\>i\<leqslant\>n,> is a symmetrical group, acting on the
    blocks of corresponding BIBD <math|B<lsub|i>.>

    \ 

    Finally, we also need to take into account the obvious combinatorial
    symmetry of the CBIBD, components that have the same parameter
    \<lambda\>. To do this, we define by <math|\<Lambda\><lsub|j> the
    \ following \ set of components:<text|>>

    <\equation*>
      \<Lambda\><lsub|j>=<around*|{|B<lsub|i><around*|\||<large|>\<lambda\><lsub|i>|\<nobracket\>>=j|}>
      for j \<in\><around*|{|\<lambda\><lsub|i><around*|\||1\<leqslant\>i\<leqslant\>n|\<nobracket\>>|}>.
    </equation*>

    direct product of the groups:

    <\equation*>
      <text|<em|<math|S<around*|(|\<Lambda\>|)>=>>>S<lsub|<around*|\||\<Lambda\>|\<nobracket\>><lsub|1><around*|\|||\<nobracket\>>>\<times\>S<lsub|<around*|\||\<Lambda\>|\<nobracket\>><lsub|2><around*|\|||\<nobracket\>>>\<times\>\<ldots\>\<times\>S<lsub|<around*|\||\<Lambda\>|\<nobracket\>><lsub|j<lsub|max>><around*|\|||\<nobracket\>>>,
    </equation*>

    where <math|S<lsub|<around*|\||\<Lambda\>|\<nobracket\>><lsub|j><around*|\|||\<nobracket\>>>=S<lsub|<around*|\||\<Lambda\>|\<nobracket\>><lsub|j><around*|\|||\<nobracket\>>><around*|(|\<Lambda\><lsub|j>|)>,><math|1\<leqslant\>j\<leqslant\>j<lsub|max>=<around*|\||<around*|{|\<Lambda\><lsub|j>|}>|\|>>,
    is a symmetrical group, acting on the corresponding components of CBIBD.

    \;

    Thus, the full set of isomorphisms of the CBIBD will be defined as a
    direct product:

    <\equation*>
      \<Gamma\><around*|(|D|)>=S<lsub|v><around*|(|V|)>\<times\><em|<text|S><around*|(|\<bbb-B\>|)>\<times\>S<around*|(|\<Lambda\>|)>>
    </equation*>

    of the groups <math|S<lsub|v><around*|(|V|)>,<em|<text|S><around*|(|\<bbb-B\>|)>,S<around*|(|\<Lambda\>|)>>>
    acting on the corresponding CBIBD elements.

    \;

    The following example of the CBIBD with parameters (6, 3, {4, 2, 2})
    illustrates the importance of considering
    <math|S<around*|(|\<Lambda\>|)>>. For that design <math|D<lsub|1>> the
    matrix

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<bbb-M\><around*|(|D<lsub|1>|)>=<around*|[|<with|font-base-size|12|<below|<with|font-base-size|14|<above||<with|font-base-size|14|<with|font-base-size|16|3\<ldots\>3<space|1em>2\<ldots\>2<space|1em>1.
      \ . . 1><rsub|>>>>|<with|font-base-size|18|<with|font-shape|italic|<em|B<rsub|1><space|1em>><space|1em><em|B<rsub|2>>
      <space|2em><em|B<rsub|<em|<em|3>>> >>>>>|]>>|<cell|>|<cell|>>>>
    </eqnarray*>

    has the form:

    <with|font-family|tt|<\with|font|roman>
      <space|9em>33333333333333333333 2222222222 1111111111

      <space|9em>11111111110000000000 1111100000 1111100000

      <space|9em>11110000001111110000 1100011100 1100011100

      <space|9em>11001100001100001111 <underline|0011011010>
      <underline|1010010011><space|1em>

      <space|9em>00101011000011101110 1000110011 0101001011

      <space|9em>00010010111010011101 0010101101 0001110110

      <space|9em>00000101110101110011 0101000111 0010101101
    </with>>

    \;

    It is easy to see that by changing the places of <math|B<lsub|2 >> and
    <math|B<lsub|3 >> we will get combined BIBD <math|D<lsub|2 >> which is
    isomorphic to <math|D<lsub|1>, \ but its representation matrix >

    <\equation*>
      <tabular|<tformat|<table|<row|<cell|\<bbb-M\><around*|(|D<lsub|2>|)>=<around*|[|<with|font-base-size|12|<below|<with|font-base-size|14|<above||<with|font-base-size|14|<with|font-base-size|16|3\<ldots\>3<space|1em>2\<ldots\>2<space|1em>1.
      \ . . 1><rsub|>>>>|<with|font-base-size|18|<em|B<rsub|1><space|1em><space|1em>><em|B<rsub|3>><space|1em><space|1em><em|B<rsub|<em|<em|2>>>
      >>>>|]>>>>>>
    </equation*>

    is lexicographically greater than the matrix
    <math|\<bbb-M\><around*|(|D<lsub|1>|)>:>
    <math|\<bbb-M\><around*|(|D<lsub|2>|)>
    \<succ\>\<bbb-M\><around*|(|D<lsub|1>|)>.> No other isomorphic
    transformation of <math|D<lsub|1> from >group
    <math|S><rsub|<em|v>>(<em|V>)<math|\<times\>><math|<text|S><around*|(|\<bbb-B\>|)>
    will give us this result. >

    <section|Enumeration results of some Combined BIBDs>

    For the enumeration of CBIBDs with parameters \ (<math|v,
    k,<around*|{|\<lambda\><lsub|1>,\<lambda\><lsub|2>,\<ldots\>\<lambda\><lsub|n>|}>>),
    we slightly modified the algorithm described in [3]. The most significant
    changes affected:\ 

    <\enumerate>
      <item>computation and storage of solutions for construction of
      <em|j>-th rows of matrices <math|<around*|{|B<lsub|i>|}>,1\<leqslant\><no-break>i\<leqslant\><no-break>n;>

      <item>checking the canonicity of fully and partially constructed
      matrices <math|\<bbb-M\><around*|(|D|)>,defined by
      <around*|(|2.2|)>,including >the use of
      <math|S<around*|(|\<Lambda\>|)>>.
    </enumerate>

    The program that implements this algorithm and some results of its use
    for enumeration of CBIBDs, BIBDs, t-designs, and semisymmetric graphs can
    be downloaded from [4].

    \;

    The following table contains some enumeration results of CBIBDs on
    <math|v elements for 6\<leqslant\><no-break>v\<leqslant\>13:<text|<strong|<next-line><next-line>Table
    1><strong|>>.>\ 

    \ <em|<space|2em><math|<around*|(|v,k,<around*|{|\<lambda\><lsub|1>,\<lambda\><lsub|2>,\<ldots\>|}>|)>>
    \ \ \ <space|2em><space|2em>Total #:<space|1em>Simple #: \ <space|2em>
    Run Time (sec):>

    <\bothlined>
      <\with|font-family|tt>
        <\with|font-base-size|9>
          <\with|font-family|tt>
            <space|3em>(6, 3, {2, 2}) \ \ \ \ \ \ \ \ \ <space|4em>2
            \ \ \ \ \ \ \ \ \ \ 2 \ \ \ \ \ \ \ \ \ \ \ \ 0.00

            <space|3em>(6, 3, {4, 2}) \ \ \ \ \ \ \ \ <space|4em>12
            \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ \ \ \ 0.02

            <space|3em>(6, 3, {4, 4}) \ \ \ \ \ \ \ \ <space|4em>24
            \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ \ \ \ 0.02

            <space|3em>(6, 3, {2, 2, 2}) \ \ \ \ \ \ <space|4em>5
            \ \ \ \ \ \ \ \ \ \ 5 \ \ \ \ \ \ \ \ \ \ \ \ 0.00

            <space|3em>(6, 3, {4, 2, 2})<space|7em>43 \ \ \ \ \ \ \ \ \ \ 2
            \ \ \ \ \ \ \ \ \ \ \ \ 0.02

            <space|3em>(6, 3, {4, 4, 2}) \ \ \ \ <space|4em>131
            \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ \ \ \ 0.04

            <space|3em>(6, 3, {2, 2, 2, 2}) \ \ <space|4em>10
            \ \ \ \ \ \ \ \ \ 10 \ \ \ \ \ \ \ \ \ \ \ \ 0.02

            <space|3em>(6, 3, {4, 2, 2, 2}) \ <space|4em>124
            \ \ \ \ \ \ \ \ \ \ 5 \ \ \ \ \ \ \ \ \ \ \ \ 0.05

            <space|3em>(6, 3, {4, 4, 2, 2}) \ <space|4em>604
            \ \ \ \ \ \ \ \ \ \ 2 \ \ \ \ \ \ \ \ \ \ \ \ 0.20

            <space|3em>(6, 3, {4, 4, 4, 2}) <space|4em>1722
            \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ \ \ \ 1.05

            <space|3em>(6, 3, {4, 4, 4, 4}) <space|4em>2490
            \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ \ \ \ 6.43

            <space|3em>(7, 3, {1, 1}) \ \ \ <space|1em>
            \ \ \ \ \ <space|3em>2 \ \ \ \ \ \ \ \ \ \ 2
            \ \ \ \ \ \ \ \ \ \ \ \ 0.00

            <space|3em>(7, 3, {2, 1}) \ \ \ \ \ <space|1em> \ \ <space|3em>16
            \ \ \ \ \ \ \ \ \ \ 3 \ \ \ \ \ \ \ \ \ \ \ \ 0.01

            <space|3em>(7, 3, {1, 1, 1}) \ \ \ \ <space|1em> \ <space|3em>5
            \ \ \ \ \ \ \ \ \ \ 5 \ \ \ \ \ \ \ \ \ \ \ \ 0.00

            <space|3em>(7, 3, {2, 1, 1}) \ \ \ \ \ <space|4em>88
            \ \ \ \ \ \ \ \ \ 20 \ \ \ \ \ \ \ \ \ \ \ \ 0.02

            <space|3em>(7, 3, {2, 2, 1}) \ \ \ \ <space|4em>607
            \ \ \ \ \ \ \ \ \ 31 \ \ \ \ \ \ \ \ \ \ \ \ 0.12

            <space|3em>(7, 3, {2, 2, 2}) \ \ \ <space|4em>2319
            \ \ \ \ \ \ \ \ \ 33 \ \ \ \ \ \ \ \ \ \ \ \ 1.38

            <space|3em>(8, 4, {3, 3}) \ \ \ <space|6em>601
            \ \ \ \ \ \ \ \ 601 \ \ \ \ \ \ \ \ \ \ \ \ 0.18

            <space|3em>(8, 4, {6, 3}) \ \ <space|4em>10648675
            \ \ \ \ \ 776249 \ \ \ \ \ \ \ \ \ 5:18.60

            <space|3em>(8, 4, {3, 3, 3}) <space|4em>1044344 \ \ \ \ 1044344
            \ \ \ \ \ \ \ \ \ 3:23.01

            <space|3em>(9, 3, {1, 1}) \ \ \ \ \ \ \ \ \ <space|4em>3
            \ \ \ \ \ \ \ \ \ \ 3 \ \ \ \ \ \ \ \ \ \ \ \ 0.03

            <space|3em>(9, 3, {2, 1}) \ \ \ \ \ \ <space|4em>8039
            \ \ \ \ \ \ \ 3514 \ \ \ \ \ \ \ \ \ \ \ \ 0.57

            <space|3em>(9, 3, {2, 2}) \ \ <space|4em>16407107 \ \ \ \ 4068169
            \ \ \ \ \ \ \ \ 14:34.27

            <space|3em>(9, 3, {3, 1}) \ \ <space|4em>15574414
            \ \ \ \ \ 208362 \ \ \ \ \ \ \ \ \ 8:09.02

            <space|3em>(9, 3, {1, 1, 1}) \ \ \ \ <space|3em> <space|1em>61
            \ \ \ \ \ \ \ \ \ 61 \ \ \ \ \ \ \ \ \ \ \ \ 0.16

            <space|3em>(9, 3, {2, 1, 1}) <space|4em>2768644 \ \ \ \ 1156282
            \ \ \ \ \ \ \ \ \ 2:51.29

            <space|3em>(9, 3, {1, 1, 1, 1})<space|4em>22727 \ \ \ \ \ \ 22727
            \ \ \ \ \ \ \ \ \ \ \ 26.06

            <space|3em>(9, 3, {2, 1, 1, 1})<space|2em>672239828 \ \ 300916861
            \ <space|2em>31:13:43.60
          </with>

          <space|3em>(10,4, {2, 2})<space|1em> \ \ \ <space|3em>
          <space|1em>7613 \ \ \ \ \ \ \ 7613 \ \ \ \ \ \ \ \ <space|1em>
          \ 4.79

          <space|3em>(10,4, {2, 2, 2})<space|1em> \ \ \ <space|3em>7613
          \ \ \ \ \ \ \ 7613 \ \ \ \ \ \ <space|3em>4.79

          <space|3em>(10,4, {2, 2})<space|1em> \ \ \ <space|3em>
          <space|1em>7613 \ \ \ \ \ \ \ 7613 \ \ \ \ \ \ \ \ <space|1em>
          \ 4.79

          <space|3em>(11,5, {2, 2})<space|1em> \ \ \ <space|3em>
          <space|1em>7613 \ \ \ \ \ \ \ 7613 \ \ \ \ \ \ \ \ <space|1em>
          \ 4.79

          <space|3em>(11,5, {2, 2, 2})<space|1em> \ \ \ <space|3em>7613
          \ \ \ \ \ \ \ 7613 \ \ \ \ \ \ \ \ <space|1em> \ 4.79

          <space|3em>(13,3, {1, 1})<space|1em> \ \ \ <space|3em>
          <space|1em>7613 \ \ \ \ \ \ \ 7613 \ \ \ \ \ \ \ \ <space|1em>
          \ 4.79

          <space|3em>(13,4, {1, 1})<space|1em> \ \ \ <space|3em>
          <space|1em>7613 \ \ \ \ \ \ \ 7613 \ \ \ \ \ \ \ \ <space|1em>
          \ 4.79
        </with>

        \;
      </with>
    </bothlined>

    In the <em|<strong|>> column \P<em|Total> #\Q of this table, we list the
    number of pairwise non-isomorphic combined BIBDs, and the column
    \ <em|<em|\P<em|Simple> #\Q<em|>>> contains the number of non-isomorphic
    CBIBDs, all of which have no replicated blocks.

    \;

    For any set of CBIBD parameters from <strong|Table 1> in [4] one can find
    detailed information about the automorphism groups of constructed sets of
    designs in the following form:

    \;

    <strong|Table 2.>

    <\with|font-base-size|9>
      <\with|font-family|tt>
        CBIBD(7, 3, {2, 1, 1})

        \;

        \ \ \ \ \ \ \ \ \|Aut(<math|D>)\| \ \ \ \ \ \ \ \ \ \ \ Nd:
        \ \ \ \ \ \ \ \ \ \ Ns: \ \ Ndt: \ \ Nst:

        \ \ \ \ _____________________________________________________________

        \ \ \ \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ \ \ \ \ \ 19
        \ \ \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 1*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 7
        \ \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 2 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 5
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 2*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 3
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 3 \ \ \ \ \ \ \ \ \ \ \ \ \ \ 17
        \ \ \ \ \ \ \ \ \ \ \ 7 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 3*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 7
        \ \ \ \ \ \ \ \ \ \ \ 2 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 4 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 4
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 4*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 6
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 6 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 6*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 1
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ \ 8*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 1
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ 12 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 2
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ 12*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 4
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ 21 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 1 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ 21*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 3
        \ \ \ \ \ \ \ \ \ \ \ 2 \ \ \ \ \ 3 \ \ \ \ \ 2

        \ \ \ \ \ \ \ \ \ \ \ \ 24 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ \ 24*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 4
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 0 \ \ \ \ \ 0

        \ \ \ \ \ \ \ \ \ \ \ 168*2 \ \ \ \ \ \ \ \ \ \ \ \ \ 1
        \ \ \ \ \ \ \ \ \ \ \ 0 \ \ \ \ \ 1 \ \ \ \ \ 0

        \ \ \ \ _____________________________________________________________

        \ \ \ \ \ \ \ \ Total: \ \ \ \ \ \ \ \ \ \ \ \ \ \ 88
        \ \ \ \ \ \ \ \ \ \ 20 \ \ \ \ \ 5 \ \ \ \ \ 2
      </with>

      <\with|font-family|tt>
        \;

        \ \ \ \ \ \ \ \ 88 matrices were constructed in \ 0.02 sec,

        \ \ \ \ \ \ \ \ \ 5 of them are transitive on the rows.

        \ \ \ \ \ \ \ \ 20 matrices have no replicated blocks,

        \ \ \ \ \ \ \ \ \ 2 of them are transitive on the rows.

        \ \ \ \ \ \ \ 187 matrices were fully constructed.

        \;
      </with>
    </with>

    If information about \|Aut(<math|D>)\| is represented as a number
    <math|N>, this means that <math|N> is the order of the group acting on
    the elements of the design, and there are no nontrivial automorphisms
    acting on the components of the corresponding CBIBDs. If it is
    represented as <em|N*M><math|>, then<math| N> and<math| M> are the orders
    of automorphism groups acting respectively on the elements and components
    of the corresponding CBIBDs.\ 
  </with>

  \;

  <section|Master-BIBD of Combined BIBD>

  As we already know, any CBIBD with parameters \ (<math|v,
  k,<around*|{|\<lambda\><lsub|1>,\<lambda\><lsub|2>,\<ldots\>\<lambda\><lsub|n>|}>>)
  is also a \Pregular\Q BIBD with parameters \ (<math|v,
  k,\<lambda\><lsub|1>+\<lambda\><lsub|2>+\<ldots\>+\<lambda\><lsub|n>>). In
  other words, if <math|D> is a CBIBD represented by the matrix
  <math|\<bbb-M\><around*|(|D|)> defined by <around*|(|3.1|)>,then >

  <\equation>
    \<bbb-B\><around*|(|D|)>=<around*|[|B<lsub|1>,B<lsub|2>,\<ldots\>,B<lsub|n>|]>
  </equation>

  is an incidence matrix of some BIBD. In what follows, we will call such
  BIBD <math|\<bbb-B\><around*|(|D|)> ><strong|a master BIBD> of <math|D. >We
  think it would be very interesting to investigate how often two or more
  non-isomorphic CBIBDs have isomorphic master BIBDs. Or a similar but
  somewhat differently worded question: how often does BIBD have different
  (non-isomorphic) representations by CDIBDs?

  <\with|par-first|0tab>
    To solve this problem, we have used the following approach.
  </with>

  \;

  <\with|par-first|0tab>
    <with|par-first|0tab|We have created a database <math|\<frak-D\> >for
    storing the canonical master-BIBDs. At the beginning of the enumeration
    process of CBIBDs with parameters \ (<math|v,
    k,<around*|{|\<lambda\><lsub|1>,\<lambda\><lsub|2>,\<ldots\>\<lambda\><lsub|n>|}>>,
    this database <math|\<frak-D\> >is empty. When the program builds the
    next canonical CBIBD <math|D,we use (5.1) to create the corresponding
    <line-break>master<no-break><no-break><text|-<no-break>>BIBD
    \<bbb-B\><around*|(|D|)><infix-and>define its canonical representation
    \<bbb-C\><around*|(|D|)>=Canon<around*|(|\<bbb-B\><around*|(|D|)>|)>.
    \ >Using \|Aut(<math|\<bbb-C\><around*|(|D|)>>)\| and matrix
    ><math|\<bbb-C\><around*|(|D|)>>, respectively, as the primary and
    secondary keys of the database <math|\<frak-D\>>, we are trying to find
    <math|\<bbb-C\><around*|(|D|)> in ><math|\<frak-D\>.
    <with|font-shape|italic|<text|<em|If such a search was not successful, we
    add the newly constructed canonical master-BIBD to our database
    <math|\<frak-D\> >with the counter equal to 1. Otherwise, we just
    increase the counter of the corresponding master-BIBD.>>> <text|>>
  </with>

  \;

  <\with|par-first|0tab>
    The results of a master-BIBD search for a given set of combined BIBD
    parameters are presented in tables similar to the following one:
  </with>

  <strong|<next-line>Table 3.><strong|>\ 

  <with|font-family|tt|<\with|font-base-size|9>
    CBIBD(7, 3, {2, 1, 1})

    \;

    \ \|Aut(M)\|: <em|Masters: \ CBIBDs: \ \ \ \ \ \ \ Distribution:>

    ________________________________________________________

    \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ 6 \ \ \ 6

    \ \ \ \ \ \ 2 \ \ \ \ \ \ \ \ \ 6 \ \ \ \ \ \ 17 \ \ \ 2*2 + 3*3 + 4

    \ \ \ \ \ \ 3 \ \ \ \ \ \ \ \ \ 5 \ \ \ \ \ \ 16 \ \ \ 1 + 3 + 3*4

    \ \ \ \ \ \ 4 \ \ \ \ \ \ \ \ \ 3 \ \ \ \ \ \ \ 8 \ \ \ 2*2 + 4

    \ \ \ \ \ \ 6 \ \ \ \ \ \ \ \ \ 4 \ \ \ \ \ \ 11 \ \ \ 1 + 2*3 + 4

    \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ 3 \ \ \ 3

    \ \ \ \ \ 12 \ \ \ \ \ \ \ \ \ 3 \ \ \ \ \ \ \ 5 \ \ \ 1 + 2*2

    \ \ \ \ \ 16 \ \ \ \ \ \ \ \ \ 2 \ \ \ \ \ \ \ 5 \ \ \ 2 + 3

    \ \ \ \ \ 21 \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ 2 \ \ \ 2

    \ \ \ \ \ 24 \ \ \ \ \ \ \ \ \ 4 \ \ \ \ \ \ 10 \ \ \ 1 + 2*2 + 5

    \ \ \ \ \ 42 \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ 2 \ \ \ 2

    \ \ \ \ \ 48 \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ 2 \ \ \ 2

    \ \ \ \ 168 \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ 1 \ \ \ 1

    ________________________________________________________

    \ \ Total: \ \ \ \ \ \ \ 33 \ \ \ \ \ \ 88 \ \ \ MaxDecomp for master: 6
  </with>>

  \;

  <\with|par-first|0tab>
    In <strong|Table 3>, master-BIBDs corresponding to all existing pairwise
    non-isomorphic Combined BIBs with parameters (<math|v,
    k,<around*|{|\<lambda\><lsub|1>,\<lambda\><lsub|2>,\<lambda\><lsub|3>|}>>)
    = (7, 3, {2, 1, 1}) are grouped by the order of their automorphism
    groups. In the <strong|<em|Masters>> column, we indicate the total number
    of non-isomorphic master BIBDs with parameters (<math|v,
    k,<around*|\<nobracket\>|\<lambda\><lsub|1>+\<lambda\><lsub|2>+\<lambda\><lsub|3>|)>=>(7,
    3, 4) that have an automorphisms group of the corresponding order. In the
    <em|<strong|CBIDs<em|>>> column, we indicate the total number of combined
    BIBDs, each of which corresponds to some master with a given automorphism
    group order. Finally, in the <em|<strong|Distribution>> column, we
    indicate how the CBIBs from the previous column are distributed among
    their masters.
  </with>

  \;

  <\with|par-first|0tab>
    For instance, from the first row of <strong|Table 3>, we can see that
    there are 6 pairwise non-isomorphic Combined BIBDs that correspond to the
    same master-BIBD. From the third row of <strong|Table 3>, we can see that
    16 pairwise non-isomorphic Combined BIBDs correspond to 5 non-isomorphic
    master-BIBDs,\ 

    <space|1em>a) one of which could be obtained in 1 way;

    <space|1em>b) one of which could be obtained in 3 different ways;;

    <space|1em>c) three of which could be obtained in 4 different ways each.

    \;

    Here \Pway\Q means the construction of a master from some Combined BIBD,
    and \Pdifferent ways\Q means that corresponding CBIBDs are different
    (non-isomorphic).

    \;

    Thus, out of 88 non-isomorphic Combined BIBDs mentioned in <strong|Table
    2>, only 33 non-isomorphic master-BIBDs could be built, and one of them
    could be obtained in 6 different ways.

    \;
  </with>

  <with|font-family|tt|<\with|font-base-size|9>
    CBIBD(6, 3, {4, 4, 4, 4})

    \;

    \ \|Aut(M)\|: <em|Masters: \ CBIBDs: \ \ \ \ \ \ \ Distribution:>

    _______________________________________________________________________________________

    \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ 5 \ \ \ \ \ 105 \ \ \ 7 + 10 + 23 + 27 +
    38

    \ \ \ \ \ \ 2 \ \ \ \ \ \ \ \ 12 \ \ \ \ \ 495 \ \ \ 2 + 2*4 + 7 + 9 +
    2*12 + 13 + 67 + 94 + 111 + 160

    \ \ \ \ \ \ 3 \ \ \ \ \ \ \ \ 17 \ \ \ \ \ 597 \ \ \ 2 + 3*3 + 3*4 + 3*6
    + 9 + 31 + 43 + 47 + 51 + 101 + 274

    \ \ \ \ \ \ 4 \ \ \ \ \ \ \ \ \ 3 \ \ \ \ \ \ 75 \ \ \ 2 + 22 + 51

    \ \ \ \ \ \ 6 \ \ \ \ \ \ \ \ 11 \ \ \ \ \ 450 \ \ \ 1 + 2*2 + 3 + 10 +
    2*13 + 28 + 30 + 130 + 218

    \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ \ 2 \ \ \ \ \ 135 \ \ \ 3 + 132

    \ \ \ \ \ 12 \ \ \ \ \ \ \ \ 12 \ \ \ \ \ 174 \ \ \ 2*1 + 4*2 + 8 + 9 +
    13 + 15 + 57 + 62

    \ \ \ \ \ 24 \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ 358 \ \ \ 1 + 3 + 7 + 10 + 43
    + 45 + 120 + 129

    \ \ \ \ \ 36 \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ 18 \ \ \ 18

    \ \ \ \ \ 60 \ \ \ \ \ \ \ \ \ 4 \ \ \ \ \ \ 71 \ \ \ 1 + 2 + 13 + 55

    \ \ \ \ 720 \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ 12 \ \ \ 12

    _______________________________________________________________________________________

    \ \ Total: \ \ \ \ \ \ \ 76 \ \ \ \ 2490 \ \ \ MaxDecomp for master: 274
  </with>>

  \;

  if we will remove Let's go back to the \Pregular\Q BIBDs, If information
  about \|Aut(D)\| is represented as a number <math|N>, this means that
  <math|N> is the order of

  <em|<em|<\strong>
    <strong|Let's go back to the \Pregular\Q BIBDs>

    \;

    \;

    <em|<em|<strong|Definition 2.1: > <strong|>><strong|>>
  </strong>>>

  \;

  [1] <itemize|><cite*|<label|CITEREFColbournDinitz2007>Colbourn, Charles J.;
  Dinitz, Jeffrey H. (2007),<nbsp><hlink|<with|font-shape|italic|Handbook of
  Combinatorial Designs>|https://archive.org/details/handbookofcombin0000unse><nbsp>(2nd<nbsp>ed.),
  Boca Raton: Chapman & Hall/ CRC,<nbsp><hlink|ISBN|https://en.wikipedia.org/wiki/ISBN_(identifier)><nbsp><hlink|1-58488-506-8|https://en.wikipedia.org/wiki/Special:BookSources/1-58488-506-8>>

  [2] A.Faradºev, Constructive enumeration of combinatorial
  objects,<nbsp><with|font-shape|italic|Problémes combinatoires et théorie
  des graphes>, Orsay 1976, Colloq. int. CNRS No. 260 (1978) pp. 131\U135.

  [3] A.V.Ivanov, Constructive enumeration of incidence
  systems,<nbsp><with|font-shape|italic|Ann. Discrete Math.>, Vol. 26 (1985)
  pp. 227\U246.

  [4] A.Ivanov, Program and some results of constructive enumeration of
  incidence systems, \ https://github.com/andrei5055/Incidence-System-Enumeration

  [5] R.C.Read, Every one a winner or How to avoid isomorphism searh when
  cataloguing combinatorial configurations, <em|Ann. Disrete Math>. 2 (1978),
  107-120.

  \;

  There are 270,474,142 Nonisomorphi 2-(9; 4; 6) Designs Patri R. J. Osterg 
  \Ward (I do have a copy)

  Contrudict my results for 2-(9,4,6)

  \;

  Östergård, P.R.J., Kaski, P. Enumeration of 2-(9, 3, \<lambda\>) Designs
  and Their Resolutions.<nbsp><with|font-shape|italic|Designs, Codes and
  Cryptography><nbsp><with|font-series|bold|27,<nbsp>>131\U137 (2002).
  https://doi.org/10.1023/A:1016558720904

  \;

  <\with|par-first|0tab>
    <itemize|>
  </with>
</body>

<\initial>
  <\collection>
    <associate|bg-color|white>
    <associate|color|black>
    <associate|info-flag|minimal>
    <associate|magnification|1>
    <associate|page-medium|paper>
    <associate|page-screen-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|CITEREFColbournDinitz2007|<tuple|5.1|9>>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|1|2>>
    <associate|auto-3|<tuple|2|2>>
    <associate|auto-4|<tuple|2|4>>
    <associate|auto-5|<tuple|3|5>>
    <associate|auto-6|<tuple|4|6>>
    <associate|auto-7|<tuple|5|7>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|4tab>| <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.15fn>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Representation
      of the Combined BIBDs> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <with|par-left|<quote|4tab>|Proposition 2.2.
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.15fn>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Isomorphisms
      of Combined BIBDs> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Enumeration
      results of some Combined BIBDs> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Master-BIBD
      of Combined BIBD> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>