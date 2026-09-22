# STL2STEP

Riconversione "morbida" di mesh (STL, o STEP nato da mesh) in B-Rep analitica con OpenCascade.

Il principio e' uno solo: **ogni modifica alla geometria e' locale**, viene verificata subito (solido
ancora chiuso, facce valide, orientamento coerente, area e volume coerenti con la mesh) e se non
passa il controllo viene scartata. Quello che non si riesce a convertire resta tassellato com'era:
il file di uscita e' sempre coerente con quello di ingresso, al massimo e' meno "pulito". Non si
cuce (Sewing) e non si ricostruisce niente globalmente.

## Fasi

- **Fase A** (`-a`): unione delle facce complanari e degli spigoli collineari.
- **Fase B** (`-b`): fori circolari — pareti cilindriche chiuse a 360°, concave, con i bordi che
  sono cerchi esatti.
- **Fase C** (`-c`): raccordi, smussi, svasature, lamature, sfere d'angolo, bossi — ogni regione di
  faccette che sta su un cilindro / cono / sfera / toro viene sostituita dalla superficie analitica.
  Dove nessuna quadrica descrive la superficie, la regione viene ricostruita con una B-spline.

Le fasi si eseguono nell'ordine in cui sono scritte in riga di comando. Senza flag la sequenza e'
`A B C A` (l'ultima riunisce le facce piane spezzate dalle sostituzioni).

```
python refit.py pezzo.stl -a -b -c -a 0.005
```

Vedi il docstring in testa a `refit.py` per la guida completa (tolleranze, come leggere il report,
il passo finale degli archi).

## File

- `refit.py` — lo strumento principale.
- `verify.py` — misura lo scarto fra un'uscita e la mesh di partenza.
- `archi.py` — conta quante spezzate della mesh sono in realta' archi di cerchio.
- `tassellate.py` — stima quanta curvatura e' rimasta tassellata nell'uscita.
- `test0.stl` … `test10.stl` — mesh di prova, ordinate per numero di facce crescente (`test0` la
  piu' semplice, `test10` la piu' complessa).
- altri script (`align.py`, `cadfit.py`, `cadinfo.py`, `extra.py`, `icp.py`, `leak.py`,
  `mancanti.py`, `rbrep.py`, `segnali.py`, `strips.py`, `vstrips.py`) — strumenti di analisi e
  debug usati durante lo sviluppo.

## Rami

- `main` — documentazione e commenti in inglese.
- `ita` — documentazione e commenti in italiano (lingua originale di sviluppo).

## Licenza

MIT — vedi [LICENSE](LICENSE).
