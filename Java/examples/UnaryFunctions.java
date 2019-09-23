package examples;

import org.quantlib.GaussKronrodAdaptive;
import org.quantlib.UnaryFunctionDelegate;
import org.quantlib.RichardsonExtrapolation;
import org.quantlib.Brent;

public class UnaryFunctions {

    public static void main(String[] args) throws java.lang.InterruptedException {
    long beginTime = System.currentTimeMillis();

    System.out.println("Integration result " +
        new GaussKronrodAdaptive(1e-8).calculate(
            new UnaryFunctionDelegate() {
		        public double value(double x) { return Math.sin(x); }
		    }, 0.0, Math.PI));

    System.out.println("Brent Solver result " +
        new Brent().solve(
            new UnaryFunctionDelegate() {
                public double value(double x) { return Math.cos(x)-x; }
        }, 1e-8, 0.5, 0.0, Math.PI));
    
    System.out.println("Richardson Extrapolation, known order " +
        new RichardsonExtrapolation(
            new UnaryFunctionDelegate() {
                public double value(double x) { return Math.exp(1 + x); }
            }, 0.1, 1.0).getValue(2.0)
        );
        
    System.out.println("Richardson Extrapolation, unknown order" +
        new RichardsonExtrapolation(
            new UnaryFunctionDelegate() {
                public double value(double x) { return Math.exp(1 + x); }
            }, 0.1, 1.0).getValue(4.0, 2.0)
        );

    }
}
