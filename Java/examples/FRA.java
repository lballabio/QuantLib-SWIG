package examples;

import org.quantlib.Actual360;
import org.quantlib.Date;
import org.quantlib.DayCounter;
import org.quantlib.FlatForward;
import org.quantlib.Month;
import org.quantlib.Settings;
import org.quantlib.YieldTermStructureHandle;
import org.quantlib.ForwardRateAgreement;
import org.quantlib.Position;
import org.quantlib.IborIndex;
import org.quantlib.Euribor3M;

public class FRA {


    public static void main(String[] args) throws Exception {

        Date todaysDate = new Date(23, Month.May, 2006);
        Settings.instance().setEvaluationDate(todaysDate);

        Date startDate = new Date(23, Month.August, 2006);
        Position.Type type = Position.Type.Long;
        double strike = 0.02;
        double notional = 100.0;
        double riskFreeRate = 0.06;
        DayCounter dayCounter = new Actual360();

        // define the underlying asset and the yield/dividend/volatility curves
        YieldTermStructureHandle flatTermStructure =
            new YieldTermStructureHandle(new FlatForward(todaysDate, riskFreeRate, dayCounter));
        IborIndex euribor3m = new Euribor3M(flatTermStructure);

        ForwardRateAgreement myFra =
            new ForwardRateAgreement(euribor3m, startDate, type, strike, notional, flatTermStructure);
        System.out.println(myFra.amount());
        System.out.println(myFra.NPV());
    }
}

