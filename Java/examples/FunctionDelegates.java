package examples;

import org.quantlib.Brent;
import org.quantlib.Newton;
import org.quantlib.DoubleVector;
import org.quantlib.OdeFctDelegate;
import org.quantlib.GaussKronrodAdaptive;
import org.quantlib.UnaryFunctionDelegate;
import org.quantlib.BinaryFunctionDelegate;
import org.quantlib.RichardsonExtrapolation;
import org.quantlib.RungeKutta;

public class FunctionDelegates {

    public static void main(String[] args) {
        System.out.println("Integration result " +
            new GaussKronrodAdaptive(1e-8).calculate(
                new UnaryFunctionDelegate() {
                    public double value(double x) { return Math.sin(x); }
                }, 0.0, Math.PI
            )
        );
    
        System.out.println("Brent Solver result " +
            new Brent().solve(
                new UnaryFunctionDelegate() {
                    public double value(double x) { return Math.cos(x)-x; }
                }, 1e-8, 0.5, 0.0, Math.PI
            )
        );
        
        System.out.println("Newton Solver result " +
            new Newton().solve(
                new UnaryFunctionDelegate() {
                    public double value(double x) { return x*x-1.0; }
                },
                new UnaryFunctionDelegate() {
                    public double value(double x) { return 2*x; }
                },                
                1e-4, 0.25, 0.1
            )
        );
                
        System.out.println("Richardson Extrapolation, known order " +
            new RichardsonExtrapolation(
                new UnaryFunctionDelegate() {
                    public double value(double x) { return Math.exp(1 + x); }
                }, 0.1, 1.0).getValue(2.0)
            );
            
        System.out.println("Richardson Extrapolation, unknown order " +
            new RichardsonExtrapolation(
                new UnaryFunctionDelegate() {
                    public double value(double x) { return Math.exp(1 + x); }
                }, 0.1).getValue(4.0, 2.0)
            );
        
        System.out.println("One dimensional adaptive Runge-Kutta result " + 
            // y'=y  and y[0] = 1
            new RungeKutta().getValue(
                new BinaryFunctionDelegate() {
                    public double value(double x, double y) { return y; }
                }, 1.0, 0.0, 1.0
            )
        );
        
        DoubleVector startVal = new DoubleVector();
        startVal.add(0.0);
        startVal.add(1.0);

        System.out.println("Two dimensional adaptive Runge-Kutta result " + 
            // y_0'=y_1 & y_1'=-y_0 and y_0[0]=0 & y_1[0]=1
            new RungeKutta().getValue(
                new OdeFctDelegate() {                    
                    public DoubleVector value(double x, DoubleVector y) {
                         DoubleVector retVal = new DoubleVector();
                         retVal.add(y.get(1));
                         retVal.add(-y.get(0));
                        return retVal; 
                    }
                }, startVal, 0.0, 0.5*Math.PI
            ).get(0)
        );        
    }
}
