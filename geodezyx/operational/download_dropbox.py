#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 14:37:57 2023

@author: psakic

Based on 
https://practicaldatascience.co.uk/data-science/how-to-use-the-dropbox-api-with-python

The Dropbox app parameters are here
https://www.dropbox.com/developers


The important permission parameters are:
    sharing.write
    the other can be read only
    
reset the token after changing the permissions!

"""

import datetime as dt
#### Import the logger
import logging
import os

import pandas as pd

# import wget
from geodezyx import conv, utils

log = logging.getLogger(__name__)


def dropbox_connect(dropbox_access_token):
    """Create a connection to Dropbox."""

    import dropbox
    from dropbox.exceptions import AuthError

    try:
        dbx = dropbox.Dropbox(dropbox_access_token)
    except AuthError as e:
        print("Error connecting to Dropbox with access token: " + str(e))
    return dbx


def dropbox_list_files(dbx, path):
    """
    Parameters
    ----------
    dbx : dropbox.dropbox_client.Dropbox
        The dropbox_connect output object.
    path : str
        the path of your Dropbox folder.
        (the root folder with your name not included)
        e.g.
        dropbox_path = "/split_2021_c/"
        dropbox_path = "/TE/RINEX 3.04/" + str(year)
        dropbox_path = "/RINEX 3.04/" + str(year)

    Returns
    -------
    df : DataFrame
        a Pandas dataframe of files in a given Dropbox folder
        path in the Apps directory.
    files : list of dropbox.files.FileMetadata
        list of dropbox.files.FileMetadata.

    """

    import dropbox

    try:

        ### need  this file loop if n files > 100
        results_block = []
        res = dbx.files_list_folder(path)
        results_block.append(res)

        while res.has_more:
            res = dbx.files_list_folder_continue(res.cursor)
            results_block.append(res)

        files_block = [res.entries for res in results_block]

        files_list = []
        for files in files_block:
            for file in files:
                if isinstance(file, dropbox.files.FileMetadata):
                    metadata = {
                        "name": file.name,
                        "path_display": file.path_display,
                        "client_modified": file.client_modified,
                        "server_modified": file.server_modified,
                        "size": file.size,
                    }
                    files_list.append(metadata)

        df = pd.DataFrame.from_records(files_list)
        # return df.sort_values(by='server_modified', ascending=False)
        df.sort_values(by="name", ascending=True, inplace=True)
        df.reset_index(drop=True, inplace=True)

        return df, files

    except Exception as e:
        print("Error getting list of files from Dropbox: " + str(e))


def dropbox_wget_cmds(
    DFfiles,
    dbx,
    out_dir_dwnld_files,
    out_dir_cmd_list,
    cmd_list_suffix,
    run_get_temp_link=False,
    run_get_perma_link=True,
    teria_archive=True,
    start_date=dt.datetime(1980, 1, 1),
    end_date=dt.datetime(2099, 1, 1),
):
    """
    Generate URL and wget commands to download files from a dropbox
    Especially to manage Teria's RINEXs but no only.

    Parameters
    ----------
    DFfiles : DataFrame
        Files DataFrame outputed from dropbox_list_files.
    dbx : dropbox.dropbox_client.Dropbox
        The dropbox_connect output object.
    out_dir_dwnld_files : str
        directory for the future downloaded files.
    out_dir_cmd_list : str
        directory for the URL list/wget commands.
    cmd_list_suffix : str
        suffix for the URL list/wget commands files.
    run_get_temp_link : bool, optional
        request temporary link. The default is False.
    run_get_perma_link : bool, optional
        request permanent link. Recommeded. Overrides run_get_temp_link.
        The default is True.
    teria_archive : bool, optional
        manage time for the Teria archive. The default is True.
    start_date : datetime, optional
        start date. The default is dt.datetime(1980,1,1).
    end_date : datetime, optional
        end date. The default is dt.datetime(2099,1,1).

    Returns
    -------
    WgetList : list
        list of wget commands.
    UrlDF : DataFrame
        DataFrame of the download URLs.

    """

    LinkList, WgetList, UrlList = [], [], []

    for ifil, fil in DFfiles.iterrows():

        fil_name = fil["name"]
        fil_path = fil.path_display

        if teria_archive:
            year = int(fil.path_display.split("/")[-2])
            fil_out_name = str(year) + "_" + fil_name
        else:
            year = 9999
            fil_out_name = fil_name

        if teria_archive:
            date_fil = conv.doy2dt(year, int(fil_name[:3]))

            if (date_fil > end_date) or (date_fil < start_date):
                log.info("skip %s", fil_name)
                continue
            else:
                log.info("handle %s", fil_name)

        if run_get_temp_link:
            try:
                result = dbx.files_get_temporary_link(fil_path)
            except Exception as e:
                print("ERROR in get temp link", e)
                continue

            LinkTuple = result.link, fil_name, fil_path

        if run_get_perma_link:
            try:
                result = dbx.sharing_create_shared_link(fil.path_display)
            except Exception as e:
                print("ERROR in get perma link", e)
                continue

            LinkTuple = result.url, fil_name, fil_path

        LinkList.append(LinkTuple)

        out_path = os.path.join(out_dir_dwnld_files, fil_out_name)

        ##### Build the wget command
        wget_cmd = " ".join(
            (
                "",
                "wget",
                "--verbose",
                "-O",
                out_path,
                "-o",
                out_path + ".log",
                #                         "-b",
                "--timeout=300",
                LinkTuple[0],
            )
        )

        print(wget_cmd)
        WgetList.append(wget_cmd)

        UrlList.append((fil_out_name, LinkTuple[0]))

        print(LinkTuple)

    #### generate the 2 URL/WGET outputfiles
    utils.create_dir(out_dir_cmd_list)

    out_cmd_list = (
        out_dir_cmd_list
        + "/"
        + "wget_cmd_"
        + str(year)
        + "_"
        + cmd_list_suffix
        + ".list"
    )
    with open(out_cmd_list, "w+") as F:
        F.write("\n".join(WgetList))

    UrlDF = pd.DataFrame(UrlList)
    out_url_list = (
        out_dir_cmd_list
        + "/"
        + "url_list_"
        + str(year)
        + "_"
        + cmd_list_suffix
        + ".csv"
    )
    UrlDF.to_csv(out_url_list, header=["fout", "url"], index=False)

    return WgetList, UrlDF, out_cmd_list, out_url_list
