o2::data::Stack::TransportFcn
transport(double etamin, double etamax)
{
  return [etamin, etamax](const TParticle& p, const std::vector<TParticle>& particles) -> bool {
           auto eta = p.Eta();
	   return (eta > etamin && eta < etamax);
	 };
}
