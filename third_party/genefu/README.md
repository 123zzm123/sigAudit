# third_party/genefu —— 26a 所需的 genefu 数据对象（只读副本）

`scripts/26a_extract_public_signatures.R` 需要从 `genefu` 的数据对象里取
已发表乳腺癌签名的**基因列表**。这里只保存**实际用到的 6 个** `.rda`，
使复现不必依赖作者本机的临时目录，也不必安装整个 genefu 包。

| 文件 | md5 | 字节 |
|---|---|---|
| `pam50.robust.rda` | `006fe5f002d9e7362f9f3b6f98ec7e44` | 2797 |
| `sig.ggi.rda` | `303016a3e7e62eb489493c3796545780` | 8209 |
| `sig.gene70.rda` | `ab60678d71e1657c550ff7ab0b994f6b` | 4729 |
| `sig.gene76.rda` | `31393d512f0805f8518cd10b6ac369dc` | 6079 |
| `sig.endoPredict.rda` | `f3e13342b92aa24ae1921e4ca53c7b5e` | 558 |
| `sig.oncotypedx.rda` | `2a9d24d90307629d1202dbfde39397b1` | 695 |

## 来源

- 包：`genefu` 2.44.0（Bioconductor）
- 获取：`curl -sL https://bioconductor.org/packages/release/bioc/src/contrib/genefu_2.44.0.tar.gz | tar xz`
  然后只取 `genefu/data/` 下的上述 6 个 `.rda`。
- **注意**：审查指出 `release/.../固定旧版本.tar.gz` **不是永久可用的获取方式**
  （Bioconductor 只保留当前 release 的源码包）。因此上表的 md5 才是核验依据：
  只要 md5 一致，就说明拿到的是同一份对象。若将来该 URL 失效，
  可用 `BiocManager::install("genefu")` 的对应版本，或直接从 Bioconductor 的
  版本归档（`.../packages/3.xx/bioc/src/contrib/`）取同版本源码包。

## 用法

```bash
SIGAUDIT_GENEFU_DATA=/path/to/third_party/genefu/data Rscript scripts/26a_extract_public_signatures.R
```

未设置该环境变量时，脚本默认使用 `<项目根>/third_party/genefu/data`。
