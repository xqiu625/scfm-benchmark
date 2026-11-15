# Draft email to GEO — GSE167363 attribute correction

**To:** geo@ncbi.nlm.nih.gov
**Subject:** GSE167363 — submitter request to correct Sex/age attributes of four samples

Dear GEO curation team,

I am the submitter of record for **GSE167363** ("Dynamic changes in human single-cell
transcriptional signatures during fatal sepsis", Qiu et al., *J Leukoc Biol* 2021,
PMID 34558746). I would like to request a correction to the sample metadata of four
samples in this series.

**Issue.** The `Sex` and `age range` attributes of the two non-survivor sepsis
patients were interchanged at submission: they are currently bound to the wrong
sample titles. The sample titles themselves are correct as deposited and match the
nomenclature used in the publication; only the Sex and age-range attributes of the
four samples below need to be corrected. All other attributes (donor, disease state,
source name, tissue) are correct and should remain unchanged.

| Sample | Title | BioSample | Current (incorrect) | Correct |
|---|---|---|---|---|
| GSM5102902 | NS LS T0 | SAMN18039640 | Sex: Male, age range: 90–95 | Sex: **Female**, age range: **65–70** |
| GSM5102903 | NS LS T6 | SAMN18039639 | Sex: Male, age range: 90–95 | Sex: **Female**, age range: **65–70** |
| GSM5511351 | NS ES_T0 | SAMN20695241 | Sex: Female, age range: 65–70 | Sex: **Male**, age range: **90–95** |
| GSM5511352 | NS ES_T6 | SAMN20695242 | Sex: Female, age range: 65–70 | Sex: **Male**, age range: **90–95** |

**Evidence (from the linked publication).** In Qiu et al. 2021, Table 1 lists the
two sepsis non-survivors as a 90–95-year-old male who died within 30 days of
enrollment and a 65–70-year-old female who died on day 1. The Results text states
explicitly: "the female nonsurvivor (P50) showed... fatality within 24 h...
Therefore, P50 was designated nonsurvivor, late-stage sepsis (NS LS) while P34 was
labeled nonsurvivor, early-stage sepsis (NS ES)." The NS LS samples are therefore
the 65–70-year-old female patient, and the NS ES samples the 90–95-year-old male
patient — the reverse of the attributes currently displayed.

Please propagate the correction to the linked BioSample records listed above. I am
happy to provide the original clinical records or any further verification needed.

Thank you for your help,
Xinru Qiu
Division of Biomedical Sciences, School of Medicine, University of California, Riverside

---

## 发送注意（不进邮件正文）

- 必须用 **GEO 提交账号关联的邮箱**发送（或抄送该邮箱）；如果当年提交账号
  在 PI（Meera Nair / Adam Godzik）名下，由他们发或让他们回复确认。
- GEO 一般 1–2 周处理 submitter 更正；只改属性、不动 accession 和标题。
- 更正后记得同步：`deliverable3` / finding 文档里标注的「GEO 当前记录（错误）」
  状态改为已更正；cohort 表的 fatal-case timing 按正确绑定写
  （NS LS = 女 65–70，入组 24h 内死亡；NS ES = 男 90–95，30 天内死亡）。
