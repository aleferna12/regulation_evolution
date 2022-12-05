LOOK AT
NOTE BELOW
:
DH is
an integer
.
It did
create problems
- not
sure why
.
The problem
was that
cells all
tend to
move to
south west
when in
group.
MAYBE fixed
-
turned DH
to a
double.


To Do
list for project
-
make patches
of food
in two
ways:
1)
a patch
stays in

place longer(so that cells have time to react to that)

2)
food self
replicates like
e.
coli (food
grows adjacent
to food
)

NOTE:
DH is
an integer
.
This can
create problems
when multiplying
with a
factor in[0, 1]
as we
do now...
-
make movement
as an
evolvable response
to food
1)
x make
energy dependent
food (allocate
particles to
maintenance or movement)
-
let decision
whether mu = startmu( or whatever, maybe
evolvable)
be in
dish
-
also how
much movement
happens should
perhaps be
a function
of cell
size?
2) -
let allocation
to resources
vs movement
be evolvable, anda
function of
food
3) -
let that
function evolve
itself,
4) -
There should
be the
possibility of

having multiple(at least two)

key/locks

- TEST: if
target area
too large, it
does not
grow anymore
!
But food
IS consumed

-
Make MAX
particles in

Dish::CellsEat()

a parameter
of class cell
-
Make possibility
for
spatial gradients
IN THE
FUTURE there
could be
some spatial
gradient?
-
Make backups, especially
to locally
restart simulations
and
make pretty
videos.
This is
a big
thing, and
should include
the following
:
-
there is
an initial
"no migration" phase,
which should
be an
alternative function
when you
start from
scratch
        but
you should
not
screw up
the time
counter
-
let command
line argument
and
parameter file
handle backup(
or
ONLY command
line?)
-
write function
that outputs
complete info
about the
field
-
write function
that begins
simulation from
backup

-
date of
birth might
be buggy(is
it still
?)

NO -
Implement Ioana
's method for persistent walk... and read the paper!
-> renske's criticism:
computational load
connectivity
        big
cells to
see effect
-
introduce biased
sustained random
walk, esp for
predators
        NO
-
introduce a
life span
parameter,
that tells
if
a cell
must die(e.g.due
to oxydative
stress or whatever...)
could be
re-
zeroed everytime
a cell
reproduces... I
guess?
- [SOLVED?]
there is
something odd
:
when predate
function is
on,
prey cells
can become
very big
for
a little
while...
does not
happen when
predate function
is off
WHY??? Is this true always? TEST!!!
-
Better random
number generator
- gitHub
-

Findcelldirections()

may make
huge mistakes
across borders,
this
can screw
up symmetry
of cell
division;

see Findcelldirections2()


DONE
- BUGS:
mother did
not
get resources
halved!
-
TWO MAJOR
OPTIMISIATIONS:
1-
Remove cells
used to
go through
the whole
field to
remove a
dead cell,
        now
it starts
from the
center of
mass andgoes
further out
until
        it
removed the
whole area
of the
cell (factor
10 optimisation).
2-

intplane update(IncreaseValEverywhere andIncreaseValIfEmpty)

used to
go throught
the whole
field andcheck
for
probability
        of
increasing food
which is
a lot
of unsuccessful
calls to
random number
generator,
now we
just ask
beforehand howmany
pixels we
update (factor
2 opt.).
-
MAJOR BUG

IN Dish::CellGrowthAndDivision2(void)

AND ALL
OF VERSIONS
BEFORE
        was
: if( c->

TargetArea()

> 2*c->half_div_area) // TARGET AREA ?????????
is now
: if( c->

Area()

> 2*c->half_div_area)
-
initialise variables
for
target vector
position correctly
-
output is
flushed when
generated (anyway,
        I
write a
lot of
lines in
one go, anda
big png
file)
-
fixed some
leaks due
to allocation
of PDE
plane (which
I am
not using)
-
food influx
location:
toggle between
everywhere (independent
of cell
position andidentity
),
only where
there are
no prey
nowhere
-> in this case
all other
parameters should
be zeroed
-
eatprob is
local to
each cell
-
Check that
key-
lock length
is extendable
to longer
lenght
(andmaybe
make a
function parser
for
trnaslatinf a
key locks
into J
values?
or
at least
something scalable
)
-
temperature might
be a
bit high
- fixed,
now it
is 20
-
Neighbour size
is back
to 2
to avoid
anisotropy,
Neigh = 2
is the
Moore neigh
-
Fixed (some
cells are
dead):
everybody should
have a
neighbour,
why is
this not true
when I
print them
?
-
mutations happen
...
-
Evolution of
adhesion as
key-
lock mechanism
:
- ca.
cpp DivideCells
should return
vector daughter_sigma
which should
be updated
with the
sigmas of
the new
cells
        andreturned
to dish
.
cpp CellGrowthanddivision2
where
        the
new
contact nrg
should be
updated for everybody for this sigma

-
Random number
generator ported
to c
++11 <random>:
somewhat tested(
for
uniqueness of
function call,
for
consistency with
seed,
graphically for randomness)
-
added predator
division area
of predators
in command
line arguments
-
read command

line argument(stump in parameters

.cpp)
(useful for
working on
cluster to
redirect output
)
- check if directory for
movie exists, exit
if
it does,
        create
it if
it does
not.
-
improve colouring
of the
field:
shade the
food where
there is
a predator
underneath
-
change name
of png
output
- Fixed:
it seems
that introduction
of medium
from outer
space is
buggy:
it happens
in a
very synchronised
manner for all (or most) cells.
-
set neighbourhood
size to
1
to avoid
predators from
eating
        beyond
contact
-
give copying
a probability
of copying
medium from
outer space(to
fix the
energy of
vacuum effect
)
-
check that
cells have
the right
neighbours all
the times<
--
THIS IS
PARAMOUNT
(expecially
when they
cross CA
boundaries)
-
checked that
contacts are
correctly assigned
after Renske
remark 21/06/18
-
area below
which a
cell is
set to
die is
now a
parameter
-
debug initcontact
length (all
kinds of
things were
buggy, e.g. infinite while loop)
-
how to
properly exit
when poulation
goes exitinct
?
Calling exit(0)
is the
right way
- check if directory for
movie exists,
if yes - do not
start the
simulation,
if no -
create andthen
start simulation
-
save data
to file
in some
way
-
make parameters
for
predation user
-defined
- Cell's ID should recycled after cell death
-
movies saved
somewhere
-
predation
        neighbour
tracking implemented
-
cells have
a map
neighbors which
works like
        {id: pair(amount of contact in pixels, duration of contact)}

cell removed
when too
small (area<
5)
- this
fixes that
issue that
full field
cells are
still present
because they
cannot disappear
