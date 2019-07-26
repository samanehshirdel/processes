import java.util.Scanner;

public class beta {


        public static final double E = 2.718281828459045;
        public static final double PI = 3.141592653589793;

        public double st_gamma(double x){
            double result= Math.sqrt(2*PI/x)*Math.pow((x/E), x);
            return result;
        }

        public static double la_gamma(double x){
            double[] p = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                    771.32342877765313, -176.61502916214059, 12.507343278686905,
                    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
            int g = 7;
            if(x < 0.5) return PI / (Math.sin(PI * x)*la_gamma(1-x));

            x = x- 1;
            double a = p[0];
            double t = x+g+0.5;
            for(int i = 1; i < p.length; i++){
                a += p[i]/(x+i);
            }

            return Math.sqrt(2*PI)*Math.pow(t, x+0.5)*Math.exp(-t)*a;
        }

        public static void main(String[] args) {
            beta test = new beta();
            System.out.print("enter a : ");
            Scanner first= new Scanner(System.in);
            double avalue = first.nextDouble();
            System.out.println(" ");
            System.out.println("enter b : ");
            Scanner second= new Scanner(System.in);
            double bvalue= second.nextDouble();
            double agamma= la_gamma(avalue);
            double bgamma= la_gamma(bvalue);
            double sumgamma= la_gamma(avalue+bvalue);
            double result= agamma * bgamma / sumgamma;
            System.out.println(" B(a,b) = " + result);





        }


        public static double sqrt(double number) {
            if (number < 0)
                return Double.NaN;
            if (number == 0 || !(number < Double.POSITIVE_INFINITY))
                return number;
            // Normalize number
            long bits = Double.doubleToLongBits(number);
            int exp = (int) (bits >> 52);
            if (exp == 0) // Subnormal x.
            {
                number *= 0x40000000000000L;
                bits = Double.doubleToLongBits(number);
                exp = (int) (bits >> 52) - 54;
            }
            exp -= 1023; // Unbias exponent.
            bits = (bits & 0x000fffffffffffffL) | 0x0010000000000000L;
            if ((exp & 1) == 1) // Odd exp, double x to make it even.
                bits <<= 1;
            exp >>= 1;

            // Generate sqrt(x) bit by bit.
            bits <<= 1;
            long q = 0;
            long s = 0;
            long r = 0x0020000000000000L; // Move r right to left.
            while (r != 0) {
                long t = s + r;
                if (t <= bits) {
                    s = t + r;
                    bits -= t;
                    q += r;
                }
                bits <<= 1;
                r >>= 1;
            }

            // Use floating add to round correctly.
            if (bits != 0)
                q += q & 1;
            return Double.longBitsToDouble((q >> 1) + ((exp + 1022L) << 52));
        }


        private static final double EXP_LIMIT_H = 709.782712893384;
        private static final double EXP_LIMIT_L = -745.1332191019411;
        private static final double LN2 = 0.6931471805599453;
        private static final double LN2_H = 0.6931471803691238;
        private static final double	LN2_L = 1.9082149292705877e-10;
        private static final double INV_LN2 = 1.4426950408889634;
        private static final double TWO_28 = 0x10000000;
        private static final double P1 = 0.16666666666666602;
        private static final double	 P2 = -2.7777777777015593e-3;
        private static final double	 P3 = 6.613756321437934e-5;
        private static final double	 P4 = -1.6533902205465252e-6;
        private static final double	 P5 = 4.1381367970572385e-8;
        private static final double SQRT_1_5 = 1.224744871391589,TWO_54 = 0x40000000000000L;
        public static double abs(double d)
        {
            return (d <= 0) ? 0 - d : d;
        }
        private static double scale(double x, int n)
        {
            if (x == 0 || x == Double.NEGATIVE_INFINITY
                    || ! (x < Double.POSITIVE_INFINITY) || n == 0)
                return x;
            long bits = Double.doubleToLongBits(x);
            int exp = (int) (bits >> 52) & 0x7ff;
            if (exp == 0) // Subnormal x.
            {
                x *= TWO_54;
                exp = ((int) (Double.doubleToLongBits(x) >> 52) & 0x7ff) - 54;
            }
            exp += n;
            if (exp > 0x7fe) // Overflow.
                return Double.POSITIVE_INFINITY * x;
            if (exp > 0) // Normal.
                return Double.longBitsToDouble((bits & 0x800fffffffffffffL)
                        | ((long) exp << 52));
            if (exp <= -54)
                return 0 * x; // Underflow.
            exp += 54; // Subnormal result.
            x = Double.longBitsToDouble((bits & 0x800fffffffffffffL)
                    | ((long) exp << 52));
            return x * (1 / TWO_54);
        }

        public static double exp(double x)
        {
            if (x != x)
                return x;
            if (x > EXP_LIMIT_H)
                return Double.POSITIVE_INFINITY;
            if (x < EXP_LIMIT_L)
                return 0;

            // Argument reduction.
            double hi;
            double lo;
            int k;
            double t = abs(x);
            if (t > 0.5 * LN2)
            {
                if (t < 1.5 * LN2)
                {
                    hi = t - LN2_H;
                    lo = LN2_L;
                    k = 1;
                }
                else
                {
                    k = (int) (INV_LN2 * t + 0.5);
                    hi = t - k * LN2_H;
                    lo = k * LN2_L;
                }
                if (x < 0)
                {
                    hi = -hi;
                    lo = -lo;
                    k = -k;
                }
                x = hi - lo;
            }
            else if (t < 1 / TWO_28)
                return 1;
            else
                lo = hi = k = 0;

            // Now x is in primary range.
            t = x * x;
            double c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
            if (k == 0)
                return 1 - (x * c / (c - 2) - x);
            double y = 1 - (lo - x * c / (2 - c) - hi);
            return scale(y, k);
        }

        private static final double TWO_52 = 0x10000000000000L;
        private static final double TWO_31 = 0x80000000L;
        private static final double TWO_64 = 1.8446744073709552e19;
        private static final double INV_LN2_H = 1.4426950216293335;
        private static final double INV_LN2_L = 1.9259629911266175e-8;
        private static final double SQRT_3 = 1.7320508075688772;
        private static final double  L1 = 0.5999999999999946;
        private static final double  L2 = 0.4285714285785502;
        private static final double  L3 = 0.33333332981837743;
        private static final double  L4 = 0.272728123808534;
        private static final double  L5 = 0.23066074577556175;
        private static final double  L6 = 0.20697501780033842;
        private static final double CP_H = 0.9617967009544373;
        private static final double  CP_L = -7.028461650952758e-9;
        private static final double  DP_L = 1.350039202129749e-8;
        private static final double CP = 0.9617966939259756;
        private static final double DP_H = 0.5849624872207642;
        private static final double OVT = 8.008566259537294e-17;
        public static int round(float f)
        {
            return (int) floor(f + 0.5f);
        }
        public static double floor(double a)
        {
            double x = abs(a);
            if (! (x < TWO_52) || (long) a == a)
                return a; // No fraction bits; includes NaN and infinity.
            if (x < 1)
                return a >= 0 ? 0 * a : -1; // Worry about signed zero.
            return a < 0 ? (long) a - 1.0 : (long) a; // Cast to long truncates.
        }
        public static double pow(double x, double y)
        {
            // Special cases first.
            if (y == 0)
                return 1;
            if (y == 1)
                return x;
            if (y == -1)
                return 1 / x;
            if (x != x || y != y)
                return Double.NaN;

            // When x < 0, yisint tells if y is not an integer (0), even(1),
            // or odd (2).
            int yisint = 0;
            if (x < 0 && floor(y) == y)
                yisint = (y % 2 == 0) ? 2 : 1;
            double ax = abs(x);
            double ay = abs(y);

            // More special cases, of y.
            if (ay == Double.POSITIVE_INFINITY)
            {
                if (ax == 1)
                    return Double.NaN;
                if (ax > 1)
                    return y > 0 ? y : 0;
                return y < 0 ? -y : 0;
            }
            if (y == 2)
                return x * x;
            if (y == 0.5)
                return sqrt(x);

            // More special cases, of x.
            if (x == 0 || ax == Double.POSITIVE_INFINITY || ax == 1)
            {
                if (y < 0)
                    ax = 1 / ax;
                if (x < 0)
                {
                    if (x == -1 && yisint == 0)
                        ax = Double.NaN;
                    else if (yisint == 1)
                        ax = -ax;
                }
                return ax;
            }
            if (x < 0 && yisint == 0)
                return Double.NaN;

            double t;
            double t1;
            double t2;
            double u;
            double v;
            double w;
            if (ay > TWO_31)
            {
                if (ay > TWO_64) // Automatic over/underflow.
                    return ((ax < 1) ? y < 0 : y > 0) ? Double.POSITIVE_INFINITY : 0;
                // Over/underflow if x is not close to one.
                if (ax < 0.9999995231628418)
                    return y < 0 ? Double.POSITIVE_INFINITY : 0;
                if (ax >= 1.0000009536743164)
                    return y > 0 ? Double.POSITIVE_INFINITY : 0;
                // Now |1-x| is <= 2**-20, sufficient to compute
                // log(x) by x-x^2/2+x^3/3-x^4/4.
                t = x - 1;
                w = t * t * (0.5 - t * (1 / 3.0 - t * 0.25));
                u = INV_LN2_H * t;
                v = t * INV_LN2_L - w * INV_LN2;
                t1 = (float) (u + v);
                t2 = v - (t1 - u);
            }
            else
            {
                long bits = Double.doubleToLongBits(ax);
                int exp = (int) (bits >> 52);
                if (exp == 0) // Subnormal x.
                {
                    ax *= TWO_54;
                    bits = Double.doubleToLongBits(ax);
                    exp = (int) (bits >> 52) - 54;
                }
                exp -= 1023; // Unbias exponent.
                ax = Double.longBitsToDouble((bits & 0x000fffffffffffffL)
                        | 0x3ff0000000000000L);
                boolean k;
                if (ax < SQRT_1_5)  // |x|<sqrt(3/2).
                    k = false;
                else if (ax < SQRT_3) // |x|<sqrt(3).
                    k = true;
                else
                {
                    k = false;
                    ax *= 0.5;
                    exp++;
                }

                // Compute s = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5).
                u = ax - (k ? 1.5 : 1);
                v = 1 / (ax + (k ? 1.5 : 1));
                double s = u * v;
                double s_h = (float) s;
                double t_h = (float) (ax + (k ? 1.5 : 1));
                double t_l = ax - (t_h - (k ? 1.5 : 1));
                double s_l = v * ((u - s_h * t_h) - s_h * t_l);
                // Compute log(ax).
                double s2 = s * s;
                double r = s_l * (s_h + s) + s2 * s2
                        * (L1 + s2 * (L2 + s2 * (L3 + s2 * (L4 + s2 * (L5 + s2 * L6)))));
                s2 = s_h * s_h;
                t_h = (float) (3.0 + s2 + r);
                t_l = r - (t_h - 3.0 - s2);
                // u+v = s*(1+...).
                u = s_h * t_h;
                v = s_l * t_h + t_l * s;
                // 2/(3log2)*(s+...).
                double p_h = (float) (u + v);
                double p_l = v - (p_h - u);
                double z_h = CP_H * p_h;
                double z_l = CP_L * p_h + p_l * CP + (k ? DP_L : 0);
                // log2(ax) = (s+..)*2/(3*log2) = exp + dp_h + z_h + z_l.
                t = exp;
                t1 = (float) (z_h + z_l + (k ? DP_H : 0) + t);
                t2 = z_l - (t1 - t - (k ? DP_H : 0) - z_h);
            }

            // Split up y into y1+y2 and compute (y1+y2)*(t1+t2).
            boolean negative = x < 0 && yisint == 1;
            double y1 = (float) y;
            double p_l = (y - y1) * t1 + y * t2;
            double p_h = y1 * t1;
            double z = p_l + p_h;
            if (z >= 1024) // Detect overflow.
            {
                if (z > 1024 || p_l + OVT > z - p_h)
                    return negative ? Double.NEGATIVE_INFINITY
                            : Double.POSITIVE_INFINITY;
            }
            else if (z <= -1075) // Detect underflow.
            {
                if (z < -1075 || p_l <= z - p_h)
                    return negative ? -0.0 : 0;
            }

            // Compute 2**(p_h+p_l).
            int n = round((float) z);
            p_h -= n;
            t = (float) (p_l + p_h);
            u = t * LN2_H;
            v = (p_l - (t - p_h)) * LN2 + t * LN2_L;
            z = u + v;
            w = v - (z - u);
            t = z * z;
            t1 = z - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
            double r = (z * t1) / (t1 - 2) - (w + z * w);
            z = scale(1 - (r - z), n);
            return negative ? -z : z;
        }

        private static final double TWO_27 = 0x8000000;
        private static final double S1 = -0.16666666666666632;
        private static final double S2 = 8.33333333332249e-3;
        private static final double	 S3 = -1.984126982985795e-4;
        private static final double S4 = 2.7557313707070068e-6;
        private static final double   S5 = -2.5050760253406863e-8;
        private static final double  S6 = 1.58969099521155e-10;
        private static final double C1 = 0.0416666666666666;
        private static final double	  C2 = -1.388888888887411e-3;
        private static final double  C3 = 2.480158728947673e-5;
        private static final double  C4 = -2.7557314351390663e-7;
        private static final double  C5 = 2.087572321298175e-9;
        private static final double  C6 = -1.1359647557788195e-11;

        private static double sin(double x, double y)
        {

            if (abs(x) < 1 / TWO_27)
                return x;  // If |x| ~< 2**-27, already know answer.

            double z = x * x;
            double v = z * x;
            double r = S2 + z * (S3 + z * (S4 + z * (S5 + z * S6)));
            if (y == 0)
                return x + v * (S1 + z * r);
            return x - ((z * (0.5 * y - v * r) - y) - v * S1);
        }




        private static double cos(double x, double y)
        {
            x = abs(x);
            if (x < 1 / TWO_27)
                return 1;  // If |x| ~< 2**-27, already know answer.

            double z = x * x;
            double r = z * (C1 + z * (C2 + z * (C3 + z * (C4 + z * (C5 + z * C6)))));

            if (x < 0.3)
                return 1 - (0.5 * z - (z * r - x * y));

            double qx = (x > 0.78125) ? 0.28125 : (x * 0.25);
            return 1 - qx - ((0.5 * z - qx) - (z * r - x * y));
        }
        public static double sin(double a)
        {
            if (a == Double.NEGATIVE_INFINITY || ! (a < Double.POSITIVE_INFINITY))
                return Double.NaN;

            if (abs(a) <= PI / 4)
                return sin(a, 0);

            // Argument reduction needed.
            double[] y = new double[2];
            int n = remPiOver2(a, y);
            switch (n & 3)
            {
                case 0:
                    return sin(y[0], y[1]);
                case 1:
                    return cos(y[0], y[1]);
                case 2:
                    return -sin(y[0], y[1]);
                default:
                    return -cos(y[0], y[1]);
            }
        }


       private static final double PIO2_1 = 1.5707963267341256, PIO2_1L = 6.077100506506192e-11, PIO2_2 = 6.077100506303966e-11,
                PIO2_2L = 2.0222662487959506e-21, PIO2_3 = 2.0222662487111665e-21,  PIO2_3L = 8.4784276603689e-32;
        private static final double TWO_20 = 0x100000 ,  TWO_16 = 0x10000, TWO_49 = 0x2000000000000L;
        private static final double TWO_24 = 0x1000000 ;
        private static final int TWO_OVER_PI[] = {
                0xa2f983, 0x6e4e44, 0x1529fc, 0x2757d1, 0xf534dd, 0xc0db62,
                0x95993c, 0x439041, 0xfe5163, 0xabdebb, 0xc561b7, 0x246e3a,
                0x424dd2, 0xe00649, 0x2eea09, 0xd1921c, 0xfe1deb, 0x1cb129,
                0xa73ee8, 0x8235f5, 0x2ebb44, 0x84e99c, 0x7026b4, 0x5f7e41,
                0x3991d6, 0x398353, 0x39f49c, 0x845f8b, 0xbdf928, 0x3b1ff8,
                0x97ffde, 0x05980f, 0xef2f11, 0x8b5a0a, 0x6d1f6d, 0x367ecf,
                0x27cb09, 0xb74f46, 0x3f669e, 0x5fea2d, 0x7527ba, 0xc7ebe5,
                0xf17b3d, 0x0739f7, 0x8a5292, 0xea6bfb, 0x5fb11f, 0x8d5d08,
                0x560330, 0x46fc7b, 0x6babf0, 0xcfbc20, 0x9af436, 0x1da9e3,
                0x91615e, 0xe61b08, 0x659985, 0x5f14a0, 0x68408d, 0xffd880,
                0x4d7327, 0x310606, 0x1556ca, 0x73a8c9, 0x60e27b, 0xc08c6b,
        };

        private static final double PI_OVER_TWO[] = {
                1.570796251296997,
                7.549789415861596e-8,
                5.390302529957765e-15,
                3.282003415807913e-22,
                1.270655753080676e-29,
                1.2293330898111133e-36,
                2.7337005381646456e-44,
                2.1674168387780482e-51,
        };

       public static double max(double a, double b)
        {

            if (a != a)
                return a;

            if (a == 0 && b == 0)
                return a - -b;
            return (a > b) ? a : b;
        }
        public static int max(int a, int b)
        {
            return (a > b) ? a : b;
        }


           private static int remPiOver2(double x, double[] y)
        {
            boolean negative = x < 0;
            x = abs(x);
            double z;
            int n;
            if (x < 3 * PI / 4) // If |x| is small.
            {
                z = x - PIO2_1;
                if ((float) x != (float) (PI / 2)) // 33+53 bit pi is good enough.
                {
                    y[0] = z - PIO2_1L;
                    y[1] = z - y[0] - PIO2_1L;
                }
                else // Near pi/2, use 33+33+53 bit pi.
                {
                    z -= PIO2_2;
                    y[0] = z - PIO2_2L;
                    y[1] = z - y[0] - PIO2_2L;
                }
                n = 1;
            }
            else if (x <= TWO_20 * PI / 2) // Medium size.
            {
                n = (int) (2 / PI * x + 0.5);
                z = x - n * PIO2_1;
                double w = n * PIO2_1L; // First round good to 85 bits.
                y[0] = z - w;
                if (n >= 32 || (float) x == (float) (w))
                {
                    if (x / y[0] >= TWO_16) // Second iteration, good to 118 bits.
                    {
                        double t = z;
                        w = n * PIO2_2;
                        z = t - w;
                        w = n * PIO2_2L - (t - z - w);
                        y[0] = z - w;
                        if (x / y[0] >= TWO_49) // Third iteration, 151 bits accuracy.
                        {
                            t = z;
                            w = n * PIO2_3;
                            z = t - w;
                            w = n * PIO2_3L - (t - z - w);
                            y[0] = z - w;
                        }
                    }
                }
                y[1] = z - y[0] - w;
            }
            else
            {
                // All other (large) arguments.
                int e0 = (int) (Double.doubleToLongBits(x) >> 52) - 1046;
                z = scale(x, -e0); // e0 = ilogb(z) - 23.
                double[] tx = new double[3];
                for (int i = 0; i < 2; i++)
                {
                    tx[i] = (int) z;
                    z = (z - tx[i]) * TWO_24;
                }
                tx[2] = z;
                int nx = 2;
                while (tx[nx] == 0)
                    nx--;
                n = remPiOver2(tx, y, e0, nx);
            }
            if (negative)
            {
                y[0] = -y[0];
                y[1] = -y[1];
                return -n;
            }
            return n;
        }

        private static int remPiOver2(double[] x, double[] y, int e0, int nx)
        {
            int i;
            int ih;
            int n;
            double fw;
            double z;
            int[] iq = new int[20];
            double[] f = new double[20];
            double[] q = new double[20];
            boolean recompute = false;

            // Initialize jk, jz, jv, q0; note that 3>q0.
            int jk = 4;
            int jz = jk;
            int jv = max((e0 - 3) / 24, 0);
            int q0 = e0 - 24 * (jv + 1);

            // Set up f[0] to f[nx+jk] where f[nx+jk] = TWO_OVER_PI[jv+jk].
            int j = jv - nx;
            int m = nx + jk;
            for (i = 0; i <= m; i++, j++)
                f[i] = (j < 0) ? 0 : TWO_OVER_PI[j];

            // Compute q[0],q[1],...q[jk].
            for (i = 0; i <= jk; i++)
            {
                for (j = 0, fw = 0; j <= nx; j++)
                    fw += x[j] * f[nx + i - j];
                q[i] = fw;
            }

            do
            {
                // Distill q[] into iq[] reversingly.
                for (i = 0, j = jz, z = q[jz]; j > 0; i++, j--)
                {
                    fw = (int) (1 / TWO_24 * z);
                    iq[i] = (int) (z - TWO_24 * fw);
                    z = q[j - 1] + fw;
                }

                // Compute n.
                z = scale(z, q0);
                z -= 8 * floor(z * 0.125); // Trim off integer >= 8.
                n = (int) z;
                z -= n;
                ih = 0;
                if (q0 > 0) // Need iq[jz-1] to determine n.
                {
                    i = iq[jz - 1] >> (24 - q0);
                    n += i;
                    iq[jz - 1] -= i << (24 - q0);
                    ih = iq[jz - 1] >> (23 - q0);
                }
                else if (q0 == 0)
                    ih = iq[jz - 1] >> 23;
                else if (z >= 0.5)
                    ih = 2;

                if (ih > 0) // If q > 0.5.
                {
                    n += 1;
                    int carry = 0;
                    for (i = 0; i < jz; i++) // Compute 1-q.
                    {
                        j = iq[i];
                        if (carry == 0)
                        {
                            if (j != 0)
                            {
                                carry = 1;
                                iq[i] = 0x1000000 - j;
                            }
                        }
                        else
                            iq[i] = 0xffffff - j;
                    }
                    switch (q0)
                    {
                        case 1: // Rare case: chance is 1 in 12 for non-default.
                            iq[jz - 1] &= 0x7fffff;
                            break;
                        case 2:
                            iq[jz - 1] &= 0x3fffff;
                    }
                    if (ih == 2)
                    {
                        z = 1 - z;
                        if (carry != 0)
                            z -= scale(1, q0);
                    }
                }

                // Check if recomputation is needed.
                if (z == 0)
                {
                    j = 0;
                    for (i = jz - 1; i >= jk; i--)
                        j |= iq[i];
                    if (j == 0) // Need recomputation.
                    {
                        int k; // k = no. of terms needed.
                        for (k = 1; iq[jk - k] == 0; k++)
                            ;

                        for (i = jz + 1; i <= jz + k; i++) // Add q[jz+1] to q[jz+k].
                        {
                            f[nx + i] = TWO_OVER_PI[jv + i];
                            for (j = 0, fw = 0; j <= nx; j++)
                                fw += x[j] * f[nx + i - j];
                            q[i] = fw;
                        }
                        jz += k;
                        recompute = true;
                    }
                }
            }
            while (recompute);

            // Chop off zero terms.
            if (z == 0)
            {
                jz--;
                q0 -= 24;
                while (iq[jz] == 0)
                {
                    jz--;
                    q0 -= 24;
                }
            }
            else // Break z into 24-bit if necessary.
            {
                z = scale(z, -q0);
                if (z >= TWO_24)
                {
                    fw = (int) (1 / TWO_24 * z);
                    iq[jz] = (int) (z - TWO_24 * fw);
                    jz++;
                    q0 += 24;
                    iq[jz] = (int) fw;
                }
                else
                    iq[jz] = (int) z;
            }

            // Convert integer "bit" chunk to floating-point value.
            fw = scale(1, q0);
            for (i = jz; i >= 0; i--)
            {
                q[i] = fw * iq[i];
                fw *= 1 / TWO_24;
            }

            // Compute PI_OVER_TWO[0,...,jk]*q[jz,...,0].
            double[] fq = new double[20];
            for (i = jz; i >= 0; i--)
            {
                fw = 0;
                for (int k = 0; k <= jk && k <= jz - i; k++)
                    fw += PI_OVER_TWO[k] * q[i + k];
                fq[jz - i] = fw;
            }

            // Compress fq[] into y[].
            fw = 0;
            for (i = jz; i >= 0; i--)
                fw += fq[i];
            y[0] = (ih == 0) ? fw : -fw;
            fw = fq[0] - fw;
            for (i = 1; i <= jz; i++)
                fw += fq[i];
            y[1] = (ih == 0) ? fw : -fw;
            return n;
        }



}
