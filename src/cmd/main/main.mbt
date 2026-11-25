///|
fn main {
  let cfg = @geo.default_config()
  let trace = @geo.run_emergence(cfg)
  println("# superemergence / MoonBit 0.10 / arena.superemergence")
  println(trace.summary())
  println("# BEGIN_CSV")
  println(trace.csv())
  println("# END_CSV")
  // ---- second pipeline: Eisenstein / Betti / ^G-operator -------------------
  let lcfg = @geo.langlands_default_config()
  let ltr = @geo.run_langlands(lcfg)
  println(
    "# langlands-emergence / Eisenstein series + non-semisimple L-parameters",
  )
  println(ltr.summary)
  println("# BEGIN_CSV_LANG")
  println(ltr.csv())
  println("# END_CSV_LANG")
  // ---- third pipeline: Frontier synthesis (M^{4|32} + B_2 Adelic + Oper + AD)
  let ftr = @geo.run_frontier_emergence()
  println(ftr.summary)
  println("# BEGIN_CSV_FRONTIER")
  println(ftr.csv())
  println("# END_CSV_FRONTIER")
  // ---- fourth pipeline: Grothendieck universes x Langlands x l-adic diamonds
  let utr = @geo.run_universe_emergence()
  println(utr.summary)
  println("# BEGIN_CSV_UNIVERSE")
  println(utr.csv())
  println("# END_CSV_UNIVERSE")
  // ---- sixth pipeline: round 6, correspondence emergence (T49-T68)
  let ctr = @geo.run_correspondence_emergence()
  println(ctr.summary)
  println("# BEGIN_CSV_CORRESPONDENCE")
  println(ctr.csv())
  println("# END_CSV_CORRESPONDENCE")
  // ---- seventh pipeline: round 7, ten coupled channels (T69-T80)
  let otr = @geo.run_octave_emergence()
  println(otr.summary)
  println("# BEGIN_CSV_OCTAVE")
  println(otr.csv())
  println("# END_CSV_OCTAVE")
  // ---- fifth pipeline: frame1-Hexad, 128 agents x 640 steps, ring vs LFM k-NN
  let hx = @geo.run_hexad_comparison()
  println(
    "# hexad-emergence / 128 agents x 640 steps / ring vs LFM2.5-Embedding k-NN",
  )
  println(hx.summary())
  println("# BEGIN_CSV_HEXAD")
  println(hx.csv())
  println("# END_CSV_HEXAD")
}
