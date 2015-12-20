auto convertLiborToContinuous(
  const auto &yield, //continuously compounded yield at point t.
  const auto &t //point t
){
  return -log(1/(1+yield*t))/t;
}
auto convertLiborToBond(
  const auto &libor, //libor rate over (0, t)
  const auto &t //t
){
  return 1.0/(libor*t+1.0);
}
auto convertBondToLibor(
  const auto &bond, //bond from 0 to t
  const auto &t //t
){
    return (1.0/bond-1.0)/t;
}
