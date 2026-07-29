

# ChiralAI

**Un motor de descubrimiento en lenguaje natural para objetivos quirales biocatalíticos: desde la consulta hasta candidatos clasificados, con rutas resueltas y verificación de viabilidad, en un solo pipeline.**

Haga la consulta en inglés. Recibirá una lista clasificada de moléculas quirales con: estereoquímica R/S definida, las clases de enzimas que las producen, una ruta biosintética predicha a partir del metabolismo central, enantioselectividad verificada por BRENDA en cada paso y viabilidad metabólica en *E. coli*, donde cada afirmación es rastreable hasta una consulta a la base de datos.

---

## El problema

~El 80 % de los nuevos fármacos quirales deben venderse como enantiómeros únicos. La biocatálisis es la ruta preferida: las enzimas son inherentemente quirales y alcanzan rutinamente >99 % de ee — pero responder a "¿qué objetivo quiral es realmente construible, por qué enzima y en qué hospedador?" hoy implica unir manualmente KEGG, BRENDA, PubChem, COBRApy y RDKit, ninguno de los cuales fue diseñado para comunicarse entre sí. ChiralAI es el tejido conectivo que faltaba.

---

## Pipeline

```
Consulta en lenguaje natural
    │
    ▼
GPT-4.1                  → 5 moléculas quirales candidatas con SMILES especificados en R/S
    │
    ▼
RDKit                    → Detección de centros estereogénicos CIP; marca centros sin asignar
    │
    ▼
KEGG                     → Consulta de vías y enzimas EC
    │
    ▼
Predictor de rutas (Nivel 1) → A* ponderado sobre reacciones KEGG; búsqueda hacia atrás
                           hasta metabolitos centrales; costos de aristas transparentes
    │
    ▼
BRENDA                   → %ee por EC, deduplicado en todos los pasos de la ruta
    │
    ▼
COBRApy / iJO1366        → FBA sobre el precursor terminal en E. coli K-12
    │
    ▼
Puntuable compuesto         → ee ponderada × Tanimoto × estereoquímica × viabilidad,
                           con nivel de confianza (high / medium / low)
    │
    ▼
CSV + JSON con marca temporal → columnas planas de puntuación + procedencia anidada
```

---

## Decisiones de diseño clave

Estas son las decisiones que hacen que la pipeline sea útil en lugar de simplemente impresionante:

- **Primero las bases de datos, con etiquetas de procedencia.** Cada valor de ee lleva la etiqueta `brenda_verified` o `llm_claim`. Las afirmaciones del LLM se descuentan un 0.6× y son visualmente distinguibles en la salida. El LLM propone; las bases de datos deciden.

- **A* ponderado con una heurística de distancia química.** El predictor de rutas no es un BFS ciego: utiliza la distancia de Tanimoto de huella Morgan al metabolito central curado más cercano como heurística, reduciendo drásticamente el espacio de búsqueda de KEGG manteniéndolo demostrablemente mejor que un enfoque voraz.

- **Costos de aristas transparentes, no una caja negra.** El costo de cada paso de la ruta se descompone en `base + termodinámico (eQuilibrator ΔG) + direccionalidad (KEGG ⇌ vs →) + anulación de reversibilidad industrial`. Un usuario puede leer el desglose y disentir con cualquier término.

- **Anulación de reversibilidad industrial.** Un ΔG puro penalizaría a las KREDs, transaminasas, IREDs, BVMOs y lipasas por funcionar "en contra de la corriente" — pero en la práctica, estas clases de enzimas se invierten rutinariamente en la biocatálisis industrial mediante reciclaje de cofactores e ingeniería de sustratos. La anulación codifica ese conocimiento experto como una búsqueda por prefijo EC.

- **Composición multiplicativa del ee.** El ee de toda la ruta es `∏(ee_i / 100) × 100` a lo largo de todos los pasos. Revela la pérdida de estereoquímica compuesta que ocultan los números de ee de un solo paso.

- **Conjunto de parada de metabolitos centrales curado.** 41 nodos: ciclo de TCA completo, glucólisis, PPP, los 20 aminoácidos, KIV. La búsqueda hacia atrás termina en realidad biológica, no en una profundidad arbitraria.

- **Tres modos de salida.** `top_n` para alternativas clasificadas, `full_tree` para el DAG de ramificación completo, `shortest_plus_diverse` para un portfolio diseñable. Los investigadores hacen diferentes preguntas; el predictor responde a las tres.

- **Caché en disco en `~/.cache/chiralai/`.** Las reacciones KEGG, archivos MOL y valores de ΔG de eQuilibrator se almacenan en caché con un TTL de 30 días. Las ejecuciones repetidas son prácticamente gratuitas.

---

## Validación

El predictor de rutas se evalúa empíricamente en comparación con la biosíntesis documentada en la literatura sobre un conjunto curado de 10 objetivos quirales (bloques de construcción quirales industriales, intermediarios farmacéuticos, metabolitos secundarios vegetales). Consulte [`benchmarks/REPORT.md`](benchmarks/REPORT.md) para ver los resultados completos por objetivo y las instrucciones de reproducción.

**Ejecución más reciente (2026-05-17):**

| Métrica | Valor |
|---|---|
| Rutas encontradas | 90 % (9/10) |
| Precursor terminal coincide con la literatura | 20 % (2/10) |
| Recuento de pasos dentro de ±2 de la literatura | 70 % |

La división 90/20 es honesta: el descubrimiento de rutas es robusto, pero la ruta predicha solo llega al punto de anclaje biológicamente correcto en ocasiones. La evaluación existe específicamente para cuantificar esto y detectar regresiones: su primera ejecución impulsó una corrección que duplicó la recuperación de coincidencias terminales (10 % → 20 %) al excluir aristas del grafo mediadas por cofactores (NADH, SAM, ATP, CoA, etc.) que la búsqueda estaba explotando como centros de conectividad. Los tres modos de fallo restantes se caracterizan y nombran en REPORT.md como los siguientes objetivos de calibración.

---

## Lo que no es (todavía)

- **Retrobiosíntesis de Nivel 2 (objetivos novedosos).** La predicción de rutas actualmente cubre compuestos que KEGG ya conoce (~12k reacciones). Los objetivos fuera de KEGG requieren RetroRules SMARTS + RDKit `RunReactants` — planeado para la Sprint 2.
- **Hospedadores distintos de *E. coli*.** La viabilidad utiliza solo iJO1366. iMM904 (levadura), P. putida, etc., aún no están integrados.
- **Ingeniería de enzimas.** ChiralAI marca cuando el ee de una enzima de tipo salvaje es insuficiente y deja el diseño de variantes al usuario (Rosetta, ProteinMPNN). No intenta esto computacionalmente.

---

## El enfoque diferencial

ChiralAI no compite con las siguientes herramientas en sus funciones principales. Ejecuta la cadena que ninguna de ellas abarca: validación quiral, predicción de rutas, enantioselectividad y viabilidad del hospedador, con procedencia explícita en todo momento.

| Herramienta | Fortaleza principal | Qué omite |
|---|---|---|
| RetroBioCat | Planificación de rutas biocatalíticas | Estereoquímica / ee |
| ASKCOS | Retrosíntesis orgánica | Enzimas |
| ChemCrow | LLM + herramientas de química | Biocatálisis |
| COBRApy | FBA a escala genómica | Estereoquímica |
| BRENDA | Base de datos de enzimas | Descubrimiento |

La disciplina que ninguna de ellas aplica: distinguir claramente las especulaciones del LLM de los hechos verificados por la base de datos.

---

## Inicio rápido

```bash
pip install -r requirements.txt
cp env.example .env
# Add OPENAI_API_KEY (required), BRENDA_EMAIL + BRENDA_PASSWORD (recommended)
python3 main.py
```

Pruebe con:
- `"enantiopure amine for a beta-lactam side chain"`
- `"chiral lactone monomer for biodegradable polymers"`
- `"(R)-secondary alcohol producible in E. coli fermentation"`

Los salidas se generan como `suggestions_<timestamp>.csv` junto con un archivo JSON complementario que contiene los árboles de rutas completos y el desglose de puntuación.

---

## Arquitectura

```
ChiralAI/
├── main.py                            Orquestador — delgado
└── ChiraLLM/
    ├── query_handler.py               GPT-4.1 → 5 candidatos como JSON
    ├── chirality_checker.py           Detección de centros estereogénicos con RDKit
    ├── database_validator.py          Analizador de archivos planos REST de KEGG
    ├── route_predictor.py             Nivel 1 A* sobre KEGG; costos transparentes
    ├── brenda_client.py               BRENDA SOAP; %ee desde comentarios
    ├── feasibility_checker.py         FBA con COBRApy sobre iJO1366
    └── enantioselectivity_scorer.py   Puntuación compuesta + nivel de confianza
└── utils/file_saver.py                CSV + JSON con marca temporal
```

---

## Salida CSV

| Columna | Fuente | Notas |
|---|---|---|
| `scoring_composite_score` | Puntuable | 0–1; ee ponderada + Tanimoto + viabilidad |
| `scoring_confidence` | Puntuable | `high` / `medium` / `low` |
| `scoring_top_enzyme_ec` | BRENDA / KEGG | EC mejor clasificado |
| `scoring_top_enzyme_ee` | BRENDA / LLM | Valor de %ee |
| `scoring_top_enzyme_source` | Puntuable | `brenda_verified` o `llm_claim` |
| `scoring_stereo_confirmed` | RDKit | Verdadero solo si todos los centros tienen R/S asignado |
| `scoring_feasibility_flux` | COBRApy | mmol/gDW/h; None si no está en el modelo |
| `route_top1_step_count` | Predictor de rutas | Pasos en la ruta predicha más corta |
| `route_top1_terminal_precursor` | Predictor de rutas | ID KEGG del metabolito central alcanzado |
| `route_top1_total_cost` | Predictor de rutas | Costo total de aristas transparentes sumado |
| `route_top1_composed_ee` | Predictor de rutas + BRENDA | ee multiplicativa a lo largo de la ruta |
| `scoring_notes` | Puntuable | Procedencia y advertencias legibles por humanos |

El archivo JSON complementario conserva la estructura anidada completa: cada paso de la ruta, cada resultado de BRENDA, cada componente de costo.

---

## Limitaciones conocidas

- **Las credenciales de BRENDA son importantes.** Sin ellas, cada sugerencia recurre a un ee de `llm_claim`. El registro es gratuito; este es el paso de configuración con mayor impacto.
- **Solo iJO1366.** Los metabolitos secundarios y muchos objetivos farmacéuticos devuelven `not_in_model`.
- **Límite de cobertura de KEGG.** El Nivel 1 no puede dirigir rutas a compuestos que KEGG no conoce. El Nivel 2 cerrará esta brecha.

---

## Hoja de ruta

- [ ] **Nivel 2** — retrobiosíntesis de objetivos novedosos mediante RetroRules SMARTS + RDKit `RunReactants`
- [ ] Modelos de hospedadores distintos de *E. coli* — iMM904 (S. cerevisiae), P. putida
- [ ] BRENDA `getEngineering` — mostrar variantes conocidas de evolución dirigida por EC
- [ ] Verificación programática de consistencia nombre ↔ SMILES CIP

---

¿Preguntas o ideas? [LinkedIn](https://www.linkedin.com/in/alexeimanuel/)
