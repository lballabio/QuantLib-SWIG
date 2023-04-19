package examples;

import org.quantlib.Date;
import org.quantlib.Month;
import java.time.LocalDate;

public class Time {


    public static void main(String[] args) throws Exception {

        Date qlDate = new Date(18, Month.April, 2023);
        LocalDate localDate = LocalDate.of(2023, 4, 18);

        Date qlDateFromLocalDate = Date.of(localDate);
        LocalDate localDateFromQlDate = qlDate.toLocalDate();

        System.out.println("qlDate:              " + qlDate);
        System.out.println("qlDateFromLocalDate: " + qlDateFromLocalDate);
        System.out.println("localDate:           " + localDate);
        System.out.println("localDateFromQlDate: " + localDateFromQlDate);
    }
}

