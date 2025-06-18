#!/usr/bin/env python3

import requests
import argparse
import os
import re
from typing import Optional

from dataclasses import dataclass


@dataclass
class ReleaseInfo:
    """Release info"""
    os: str
    arch: str
    idx: int = 0


@dataclass
class Args:
    """Command-line args"""
    tag: str
    api_url: str


# --------------------------------------------------
def get_args() -> Args:
    parser = argparse.ArgumentParser(
        description="Create a release table for a GitHub release"
    )

    parser.add_argument(
        "repo", help="The github repo"
    )

    parser.add_argument(
        "tag", help="The tag of the release"
    )

    args = parser.parse_args()

    return Args(tag=args.tag, api_url=f"https://api.github.com/repos/{args.repo}/releases")


def main() -> None:
    args = get_args()
    releases = get_releases_data(args.api_url)

    if release := find_release_by_tag(releases, args.tag):
        print("tag '{}' has {} assets".format(args.tag, len(release["assets"])))
        markdown_table = generate_markdown_table(release)
        update_release_body(args.api_url, release["id"], markdown_table)
        print(f"tag'{args.tag}' updated successfully.")

    else:
        print(f"tag '{args.tag}' not found.")


def get_releases_data(url: str):
    headers = {
        "Authorization": f'token {os.getenv("GITHUB_TOKEN")}',
        "Accept": "application/vnd.github.v3+json",
    }
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    return response.json()


def find_release_by_tag(releases, tag):
    for release in releases:
        print(release)
        if release["tag_name"] == tag:
            return release
    return None


ARCH_MAP = {
    "macos-universal": ReleaseInfo("MacOS", "Universal"),
    "x86_64-apple-darwin": ReleaseInfo("MacOS", "Intel 64 bit"),
    "aarch64-apple-darwin": ReleaseInfo("MacOS", "ARM 64 bit (Apple Silicon)"),
    "x86_64-unknown-linux-gnu": ReleaseInfo("Linux", "Intel/AMD 64 bit"),
    "aarch64-unknown-linux-gnu": ReleaseInfo("Linux", "ARM 64 bit"),
    "x86_64-pc-windows-msvc": ReleaseInfo("Windows", "Intel/AMD 64 bit"),
}

# redefine arch map to give the objects
# indices in the order they are defined
ARCH_MAP = {
    k: ReleaseInfo(v.os, v.arch, idx=i)
    for i, (k, v) in enumerate(ARCH_MAP.items())
}


def extract_os_arch_from_filename(filename) -> Optional[ReleaseInfo]:
    pattern = re.compile(r"^(.+)(?:\.tar\.gz|\.tgz|\.zip)$")
    if match := pattern.search(filename):
        _, stem = match.group(1).split("-", 1)
        print(f"stem: {stem}")
        if stem in ARCH_MAP:
            return ARCH_MAP[stem]
        else:
            print(f"unexpected stem: {stem}")
    else:
        print(f"no match: {filename}")


# --------------------------------------------------
def generate_markdown_table(release) -> str:
    table = "### Release Assets\n"
    table += "| OS | Architecture | Link |\n"
    table += "|---------|----------|-------------|\n"

    artifacts = [
        (info, asset)
        for asset in release["assets"]
        if (info := extract_os_arch_from_filename(asset["name"]))
    ]

    assert (len(artifacts) == len(ARCH_MAP))

    artifacts.sort(key=lambda x: x[0].idx)

    for (info, asset) in artifacts:
        print(">>> asset {}".format(asset["name"]))
        download_url = asset["browser_download_url"]
        table += f"| {info.os}  | {info.arch}  | [Download]({download_url}) |\n"

    bin, _ = artifacts[0][1]["name"].split("-", 1)
    table += (
        f"\nOn MacOS, you'll see the following message when you try to run `{bin}`: \n"
        "\n```\n"
        f"\"{bin}\" cannot be opened because the developer cannot be verified."
        "\n```\n"
        "\nTo fix this, you can run the following command:\n"
        "\n```\n"
        "sudo xattr -dr com.apple.quarantine <path to binary>\n"
        "```\n"
    )

    return table


# --------------------------------------------------
def update_release_body(api_url, release_id, new_body):
    headers = {
        "Authorization": f'token {os.getenv("GITHUB_TOKEN")}',
        "Accept": "application/vnd.github.v3+json",
    }
    update_url = f"{api_url}/{release_id}"
    data = {"body": new_body}

    response = requests.patch(update_url, headers=headers, json=data)
    response.raise_for_status()


# --------------------------------------------------
if __name__ == "__main__":
    main()
