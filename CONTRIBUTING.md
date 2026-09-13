# 参与贡献

本仓库是 2026 年高教社杯 A 题的四套解法与交付记录。欢迎队友补充、修正与继续完善。

---

## 一、先加入组织

仓库属于组织 [`attachment-too-large`](https://github.com/attachment-too-large)：

```
https://github.com/attachment-too-large/cumcm-2026-a-attempts-and-thoughts
```

**加入方式：把用户名或邮箱发给管理员，由管理员发出邀请。**

GitHub 的组织邀请**必须逐个发出**，没有"谁点谁能进"的通用链接
（参见 [GitHub 文档](https://docs.github.com/zh/organizations/managing-membership-in-your-organization/inviting-users-to-join-your-organization)）。
你可以二选一：

- **给 GitHub 用户名** —— 管理员在 https://github.com/orgs/attachment-too-large/people
  → **Invite member** → 输入用户名 → 发送
- **给邮箱** —— 同上，输入邮箱即可。
  注意：该邮箱必须是**你 GitHub 账号上已验证的邮箱**，否则收不到也接受不了邀请

邀请发出后你会收到一封邮件，**点里面的链接**才能加入。邀请 7 天后过期，可让管理员重发。

**加入组织之后，如需直接推送权限**，管理员再加为仓库协作者（网页或 API）：

```
PUT /repos/attachment-too-large/cumcm-2026-a-attempts-and-thoughts/collaborators/<你的用户名>
Content-Type: application/json

{"permission": "write"}
```

> 说明：组织成员在 People 页面公开可见（组织可以设置成员列表可见性）。
> 本仓库此前做过匿名化处理（四套方案统一以「方案 1–4」编号），
> 加入组织意味着公开你的账号身份。如果希望保持匿名，见文末「匿名参与」。

---

## 二、本地准备

```bash
git clone https://github.com/attachment-too-large/cumcm-2026-a-attempts-and-thoughts.git
cd cumcm-2026-a-attempts-and-thoughts

# 确认提交身份（很重要：这与「贡献者」头像墙直接相关）
git config user.name  "你的名字"
git config user.email "你的GitHub邮箱"
```

**邮箱必须是绑定到你 GitHub 账号的那个**，否则你的提交不会计入仓库的贡献者列表。

---

## 三、三种提交方式

### 方式 1：小改动直接推 main

```bash
git pull
# ...修改文件...
git add -A
git commit -m "修正方案 2 的边界处理说明"
git push
```

### 方式 2：开分支 + Pull Request（改动较大时推荐）

```bash
git checkout -b fix/fangan2-boundary
# ...修改...
git commit -m "..."
git push -u origin fix/fangan2-boundary
```

然后在网页上开 Pull Request，由管理员合并。

### 方式 3：在你的提交里标记共同作者

如果一个改动是两人一起做的，在提交信息末尾加尾注：

```
git commit -m "补做方案 4 的生产路径收敛性验证

Co-authored-by: 队员B <b@example.com>
Co-authored-by: 队员C <c@example.com>"
```

GitHub 会在这条提交上显示多个头像。**同样要求那个邮箱绑定过 GitHub 账号**，否则只显示文字。

---

## 四、本仓库的几条约定

1. **文件编码 UTF-8（无 BOM），换行 LF。** 仓库已配置 `.gitattributes` 自动处理。
2. **二进制文件不参与换行转换**（pdf / png / xlsx / zip / npz 已在 `.gitattributes` 中标明）。
3. **提交信息用中文**，第一行控制在 50 字以内，需要展开就空一行再写正文。
4. **不要提交 `__pycache__`、`*.pyc`、LaTeX 编译中间件**（`.aux/.log/.out/.toc`），`.gitignore` 已覆盖。
5. **PDF 有配额约定**：本仓库目前只公开两份 PDF——最终论文与讲解教材。新增 PDF 前请先说明用途。
6. **方案编号而非姓名**：四套方案在仓库内统一称「方案 1」至「方案 4」，请沿用这套编号，不要改回代号。
7. **改数值必须说明**：任何会改变结果文件数字的改动，请在提交信息里写清改了什么、为什么，以及是否重跑了验证。

---

## 五、值得继续做的事

以下是目前已知的缺口，欢迎认领：

| 事项 | 说明 | 相关目录 |
| --- | --- | --- |
| 核对文中沿用数字 | 详解教材里若干「复跑报告」数字取自较早一次运行（`run_all.py` 未重跑），可复算核对 | `method/详解/` |
| 补全 `result2` 的全时长输出 | 目前只有方案 2 按题面输出了 0–259200 s，其余三份为 0–3 h | `comparison/方案1~3/` |
| 端面效应的量级估计 | 目前按一维径向处理，未定量估计端面影响 | `delivery/A题_药材烘干_正式交付/code/` |
| 补充二维轴对称算例 | 用于估计一维简化忽略的端面对流影响量级 | `delivery/A题_药材烘干_正式交付/code/` |

---

## 六、匿名参与

如果你希望贡献但**不想公开 GitHub 账号**，有两种办法：

1. **只提内容，不推送**：把改动发在 Issue 里，或把补丁文件给管理员，由管理员代为提交，并在提交信息里用 `Co-authored-by` 标注（若邮箱绑定账号，仍会显示头像；若不想显示，就不加这行）。
2. **在 README 里以文字署名**：在「作者与贡献者」一节按方案编号写清楚，不使用 GitHub 账号。

---

## 七、有问题找谁

仓库 Issues 页面留言即可。
