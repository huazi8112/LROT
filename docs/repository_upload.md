# Independent repository preparation

This directory contains the code, retained results/images, fixed GO/MGI inputs
and S1/S2. Put its **contents** at the Git repository root so that README.md,
data/, paper_protocol/ and the other code folders are directly visible.
No sibling delivery folders, private backup, build script or manuscript belongs
in the uploaded tree. A ZIP uploaded as a single file does not expand into code.

Before uploading, run these checks from the repository root:

```text
python tools/verify_repository_package.py
python tools/verify_inputs.py --dataset go --data-root .
```

Neither check installs packages, accesses the network, imports experimental
dependencies or runs an experiment. The first checks files, current local links,
fixed source hashes and selected attachment metadata; the second hashes GO inputs.
Original tissue inputs and CellChatDB binaries remain separately obtained under
[their source and checksum instructions](../data/README.md).

## Git transfer

Browser upload cannot transmit this complete package: it contains files above
25 MiB and more than 100 files. Use an existing Git installation or client.
No automatic uploader or credentials are supplied. Large saved inputs are part
of the allowlist and must not be silently omitted, recompressed, renamed or
replaced by LFS pointer text. The included .gitattributes disables line-ending
conversion to preserve file identities. Keep both hidden Git configuration files.

Clone the intended repository into a new local working directory. Copy only the
contents of this package into it, preserving the clone's .git directory. Verify
the copy before staging. With a Git version supporting --pathspec-from-file:

```text
git config --local core.autocrlf false
python tools/verify_repository_package.py
git add -f --pathspec-from-file=RELEASE_FILES.txt
git add -f -- MANIFEST.sha256 MANIFEST_size.txt
git diff --cached --stat
git commit -m "Prepare standalone LROT repository"
```

Inspect the staged changes before committing. A push is a separate authorized
network action. No push, tag, release, DOI or remote deletion was performed when
preparing this package. Never force-push to erase history. If the existing branch
still contains obsolete delivery wrappers, remove only the identified old wrapper
paths in that local checkout after backup; do not use a broad recursive delete.

After transfer, download/clone a fresh copy and rerun the same checks. Local
integrity does not prove remote completeness. This folder remains a technical
package, not evidence of whole-paper end-to-end reproduction or final submission
approval. Existing collective citation/version metadata are retained; no author,
software release number or archive identifier is invented.

## 中文要点

只需处理这一个目录；GO和S1已在其内部。将目录**内部全部内容**放入仓库根目录，
不要把外部的备份、修改记录、打包脚本或整个交付工作目录上传。根目录数据声明可写
`data/README.md`和`data/input_files.json`。网页不能完整传输本包的大文件，请使用
已有Git；保留`.gitattributes`、`.gitignore`及全部清单文件。先校验本地副本，核对
暂存内容，再提交；远程上传完成后重新下载校验。不要安装依赖、改实验或跳过大文件
来满足上传。文件完整不等于所有实验已经端到端复现。
