#ifndef lint
static    char sccs_id[] = "@(#)rans.c    1.1    (ucb.bach)    3/9/84";
#endif

/*
 * SUPER-DUPER 'UNI' RANDOM NUMBER GERNERATOR
 *
 * PDP-11 and VAX-11 version of SuperDuper random number
 * generator (Marsaglia, Ananthanarayanan & Paul, 1973,
 * Learmonth & Lewis, 1973, and Dudewicz, 1976).
 * Implementors:
 *        Steve Hubert     2 March 82
 *        Jim Reeds    June/July  82
 *        Mark Abrahams    5 August 83
 *
 *
 * C LANGUAGE USAGE:
 *
 *    int c, t, u, low, high;
 *    int lran(), intran();
 *    short sran();
 *    double rans();
 *    int rand();
 *
 *    lran()            returns a random signed int
 *    sran()            return a random signed short
 *    rans()            returns a random real in [0,1)
 *    intran(low, high)    returns random int in range low <= x <= high
 *    random()        initialize with physical random seed
 *    getsd(&c, &t)        put current values of two seeds in c and t
 *    setsd(c, t)        set current seeds
 *
 * FORTRAN LANGUAGE USAGE:
 *
 * (The fortran interface routines are in the file rans_.c)
 *
 *    integer*4 c, t, a, b, i, lran, intran
 *    double precision x, rans
 *
 *    call random        random starting point
 *    call getsd(c, t)    get values of seeds
 *    call setsd(c, t)    set values of seeds
 *    x = rans()        make a random real in [0,1)
 *    x = rans(dummy)        ditto
 *    i = lran()        make a random signed integer
 *    i = lran(dummy)        ditto
 *    i = intran(a, b)    make a random integer in {a,a+1,...,b}
 *
 *
 * DESCRIPTION:
 *
 * Two independent generators are maintained internally:
 * a linear congruential generator with multiplier 69069 and modulus 2^32;
 * and a Tausworthe generator based on primitive trinomial x^32 + x^17 + 1.
 * The two seeds are XOR'd together to give the random number.
 * If either seed is 0 the output consists of only the other generator.
 * If both seeds are 0 only 0's are returned.
 *
 * See the manual entry on rans or rans(3), for more details
 * and for references.
 *
 *
 * SEED INFO:
 *
 * The only requirements are that a maximal period be obtained for
 * each generator.  They are maximal provided:
 *        0 < cseed < 2^32, and is odd.
 *        0 < tseed < 2^32.
 * See Knuth, Vol 2, Section 3.2.1.2, Thm B (page 19 of 2nd edition);
 * and Tausworthe (1965), page 201, for relevant theorems.
 *
 * Gross (1978) makes certain recommendations regarding initial values of the
 * seeds; however they do not apply to this implementation of the generator.
 * Gross seems to use 2^31 as the modulus for the congruential generator;
 * however we use 2^32.
 */

static int tseed = 1073, cseed = 12345;
int lran();


/*
 * 32 bit Tausworthe (linear feedback shift register) generator:
 * Kill sign extension in the right shift.
 * Zero extension on the left shift is guaranteed by C.
 *
 * Congruential generator:
 * On PDP11/VAX series this multiplication is modulo 2^32.
 */

#define UNI \
    tseed ^= ((tseed >> 15) & 0377777); \
    tseed ^= (tseed << 17); \
    cseed *= 69069

/*
 * Return the random 32 bit pattern
 */
int lran() {
    UNI;
    return ( cseed ^ tseed);
}

/*
 * Return a uniform double precision, 0 < rans < 1
 * Repeat coding of lran() for speed
 * The hardware wants to think of the 32 bit pattern as obeying
 *         - 2^31 < lran() < 2^31
 * but the logic of the algorithm wants
 *        0 < lran() < 2^32.
 * Hence the x += 1.
 */
double rans_() {
    double x;
    UNI;
    x = ((double)(tseed^cseed)) / 4294967296.;
    if (x < 0.) x += 1.;
    return(x);
}

/*
 * Return a random integer between low and high (inclusive)
 */
int intran(low,high)
int low,high;
{
    if (high < low)
        return (0);
    else
        return (low + (int)(rans_()*(high+1-low)));
}

/*
 * random start
 */
#include    <sys/types.h>
random_() {
    register int teven;
    time_t    time();

LOOP:
    /*
     * Even if time and pid don't change between calls make the seeds
     * different.  This is the DeVogelaere propaganda bug.
     */
    cseed += time((int *)0);
    tseed += getpid();

    /*
     * Make them both odd for now; but record whether t was even.
     */
    teven = ! (tseed & 1);
    tseed |= 1;
    cseed |= 1;

    /*
     * Mix up the bits some.  This operation is invertible;
     * preserves oddness, and hence maps nothing to 0.
     * (Odd numbers have unique multiplicative inverses modulo 2^p).
     * So it won't lose randomness.
     */
    cseed = tseed + (cseed^tseed) * 314159;
    tseed = cseed + (cseed^tseed) * 1492365;

    /*
     * Restore parity of tseed - still less randomness is lost.
     * But cseed must remain odd.
     */
    tseed ^= teven;

    /*
     * Just in case the first or last steps gave 0.
     */
     if (!tseed) goto LOOP;

    /*
     * For good measure, run off a couple - it just mixes more.
     */
    lran();
    lran();
}

/*
 * Return the most significant 16 bits of the 32 bit pattern
 */
short sran() {
    return (0xffff & (lran() >> 16));
}

/*
 * Set and get the seed values
 */
setsd(c,t)
int c,t;
{
    cseed = c;
    tseed = t;
}

getsd(cp,tp)
int *cp,*tp;
{
    *cp = cseed;
    *tp = tseed;
}
