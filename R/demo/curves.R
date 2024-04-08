

suppressMessages(library(QuantLib))


bToday = QuantLib::Date(31, "August", 2022)
QuantLib::Settings_setEvaluationDate(self = QuantLib::Settings_instance(), d = bToday)
print(Settings_instance()$getEvaluationDate())

bufDates = QuantLib::DateVector()
QuantLib::DateVector_append(self = bufDates, x = Date(31, "August", 2022))
QuantLib::DateVector_append(self = bufDates, x = Date(30, "September", 2022))
QuantLib::DateVector_size(self = bufDates)

bCurve = QuantLib::DiscountCurve(dates = bufDates , discounts = c(1.0, 1.0), dayCounter = Actual360())


