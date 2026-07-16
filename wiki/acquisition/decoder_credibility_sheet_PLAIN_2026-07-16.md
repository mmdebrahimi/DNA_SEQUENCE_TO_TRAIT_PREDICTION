# dna-decode, explained without the jargon

*A plain-language companion to `decoder_credibility_sheet_2026-07-16.md`. Same facts, same numbers — no
genomics background needed. Research use only; this is not a medical device and must not guide treatment.*

---

## The problem, in one paragraph

When someone has a bacterial infection, doctors need to know which antibiotics will still kill it. The
standard way is to grow the bacteria in a lab and try the drugs on it — accurate, but it takes **days**.
Meanwhile we can read a bacterium's **entire DNA in hours**. So: can you read the DNA and *predict* which
drugs will fail, without waiting for the lab?

That's what this tool does. Same idea for HIV, tuberculosis, some fungi — and, more recently, a few
non-medical traits (a human's eye colour; whether a plant flowers early or late).

---

## How it works — and why that matters

Think of antibiotic resistance like a **lock and key**. The antibiotic is the key; a part of the bacterium is
the lock. Bacteria become resistant by **changing the lock** or by **carrying a tool that destroys the key**.

Scientists have spent decades cataloguing those specific changes: *this* mutation in *this* gene means *this*
drug stops working. Our tool does something deliberately unglamorous — **it looks for the catalogued changes
and reports what it finds.**

That's the opposite of how most modern AI works. We are **not** showing a computer thousands of examples and
letting it find its own patterns. We're applying known biology, mechanically.

**Why that's a feature, not a limitation:**

| | pattern-matching AI | this tool |
|---|---|---|
| Can it explain itself? | "the model says resistant" | "**resistant — because I found gene X at position Y**" |
| Can an expert check it? | not really | yes, line by line |
| What if it doesn't know? | guesses anyway | **says "I don't know"** |
| Does it drift over time? | yes | no — the rules are frozen and cryptographically checked on every run |

That last one matters more than it sounds. **The tool refuses to answer when the biology isn't in the
catalogue.** Most tools guess. A confident wrong answer is worse than no answer.

---

## The honest part (this is the whole point)

Anyone can claim a high accuracy number. Here's the one we *don't* hide:

> ### The photocopy problem
>
> We tested our tuberculosis predictions on **2,845 real TB samples**. We got them right **92% of the time**.
> That number is misleading, and we say so.
>
> Bacteria reproduce by **cloning themselves**. So those 2,845 samples aren't 2,845 independent tests — they
> collapse into roughly **67 genuinely different family trees**. It's like testing a spell-checker on 2,845
> documents when 2,700 of them are photocopies of the same 67 originals. You'd report "99% accurate!" while
> having really only tested 67 things.
>
> When we count **once per family** instead of once per photocopy, our score drops from **92% to 44%**.
>
> **We publish the 44%.** That's our headline number.

We do the same everywhere: we correct for it, we show the uncertainty range, and we let the ugly number lead.

**And we publish our failures — including the ones that turned into wins.**

We tried the fashionable shortcut first: take a big pre-trained AI "foundation model" **off the shelf**, don't
train it on anything, just ask it. That failed on **five different problems**, every time. In one case it did
**worse than random guessing** at a task where our boring rulebook scored 93%. We wrote up exactly why instead
of quietly dropping it.

**But that's not the whole story, and the rest of it matters more.**

The rulebook has a real blind spot: some samples are resistant with **no catalogued mutation at all**. The
rulebook simply cannot see them. So we did the thing the shortcut skipped — we **properly trained** a model on
HIV, the one problem where we have thousands of real lab measurements to learn from.

**It worked.** It catches the cases the rulebook structurally misses, and it holds up when tested on studies
it never saw. (Concretely: **0.81** on that held-out test — where the off-the-shelf AI managed 0.449, barely
better than a coin flip.)

**So what actually ships is a hybrid**, not a rulebook purist:

> **The rulebook leads** — explainable, checkable, abstains when unsure.
> **The learned model rides alongside it**, covering the rulebook's blind spot.

We also found **where that stops working**, which is the genuinely interesting bit:

- **It works for HIV** — a fast-mutating virus where the same resistance mutations keep arising independently.
  There's a real pattern to learn.
- **It fails for tuberculosis** — a slow, clonal bacterium. There, the model *looked* like it worked (0.66),
  but it was really just **recognising family trees**. Hide the families and it dropped to **0.51 — a coin
  flip**. We caught that, and did not ship it.

That last one is the whole method in miniature: **the honest test is the one that hides what the model could
memorise.** Ordinary testing said "ship it." The honest test said "this is a lineage-memoriser." We believe the
HIV result *precisely because* it survived the same kind of test that killed the TB one.

If someone shows you a genomics tool with one big impressive accuracy number and no caveats, that's the thing
to be suspicious of — not this.

---

## What's actually been built

**27 different prediction tasks** ("can you predict drug X in organism Y?"). Of those:

- **10 are properly scored** — tested on samples from *different labs in different countries* than the ones
  the tool was built from. (That's the meaningful test. Testing on data from the same source you learned from
  is like marking your own homework.)
- **2 deliberately refuse to answer** — the biology isn't decodable from genes alone, so the tool abstains.
- **11 can't be scored at all** — not because the tool fails, but because **nobody has published the lab
  results needed to check it**. Which brings us to the ask.

A couple of specific results, translated:

- **HIV, efavirenz: 0.962.** Meaning: hand it one drug-resistant virus and one susceptible one, and it puts
  them in the right order **96% of the time**. This one is special because we checked it against **actual
  wet-lab measurements** — real experiments, not another computer's opinion.
- **Tuberculosis: 44%** sensitivity after the photocopy correction (see above).

It's free, open-source, and anyone can install it: `pip install dna-decode`.

---

## What we're stuck on — and what we'd need from you

Here's the honest bottleneck, and it isn't cleverness.

**We have the methods. We have the DNA. What we don't have is the answer key.**

To check whether a prediction is right, you need two things *for the same sample*: its **DNA sequence**, and
**what the lab actually measured** when they tested the drugs on it. Those pairs are surprisingly rare in
public. Loads of DNA out there. Loads of lab results out there. Very few places where you can get **both,
for the same bacterium**.

**If you have that pairing, here's exactly what we'd do:**

1. Run our **already-finished** tool on it. We would **not** retrain or tune anything — which is the point.
   It's a genuinely fair test precisely because the tool can't have been fitted to your data.
2. Report how it did, **including the failures**, with the photocopy correction and honest uncertainty ranges.
3. Send you the analysis.

**What we'd need:** the DNA sequence, the lab result, and which method the lab used. **That's it.**
No patient names, no medical records, no personal information of any kind.

**What we wouldn't do:** share your data with anyone, or publish anything without your say-so. Happy to work
under whatever formal data-sharing agreement your institution requires.

---

## The short version

We built a DNA reader that predicts drug resistance by looking up known biology instead of guessing from
patterns — so it can show its work, and it says "I don't know" rather than bluffing. We report the numbers
that make us look worse, and we publish the experiments that failed. What we're missing isn't a better
algorithm; it's **DNA paired with real lab results** to check ourselves against.

*Every number here is generated automatically from the project's own records — the full technical version,
with sources, is the companion document.*
